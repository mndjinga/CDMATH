import MEDLoader as ml

def read_typ2(fichier, nom_med):
    with open(fichier, "r") as fic: lines = fic.readlines()
    nb_som = int(lines[1].strip())
    nb_cel = int(lines[3 + nb_som].strip())

    som = [map(float, line.split()) for line in lines[2:2 + nb_som]]

    # les connectivites ne partent pas de 0...
    tmp_cel = [map(int, line.split()[1:]) for line in lines[4 + nb_som:]]
    cel = [[p - 1 for p in c] for c in tmp_cel]

    mesh = ml.MEDCouplingUMesh("mesh", 2)
    mesh.allocateCells(len(cel))
    for p in cel: mesh.insertNextCell(ml.NORM_POLYGON, len(p), p)
    mesh.finishInsertingCells()

    pts = []
    for p in som: pts.extend(p)
    co = ml.DataArrayDouble(pts, len(pts) / 2, 2)
    mesh.setCoords(co)

    mf, d, di, r, ri = mesh.buildDescendingConnectivity()
    mf.setName("mesh")
    mm = ml.MEDFileUMesh.New()
    mm.setMeshAtLevel(0, mesh)
    mm.setMeshAtLevel(-1, mf)

    g = []
    nb_vois = ri.deltaShiftIndex()
    for i in range(mf.getNumberOfCells()):
        if nb_vois[i] == 1: g.append(i)
    grp = ml.DataArrayInt.New(g)
    grp.setName("boundary")
    mm.addGroup(-1, grp)

    mm.write(nom_med, 2)

def read_typ3(fichier, nom_med):
    with open(fichier, "r") as fic: lines = fic.readlines()
    nb_som = int(lines[9].strip())
    nb_cel = int(lines[11].strip())
    nb_fac = int(lines[13].strip())

    som = [map(float, line.split()) for line in lines[17:17 + nb_som]]


    # les connectivites ne partent pas de 0...
    tmp_cel = [map(int, line.split()[1:]) for line in lines[18 + nb_som:18 + nb_som + nb_cel]]
    cel = [[p - 1 for p in c] for c in tmp_cel]

    tmp_fac = [map(int, line.split()[1:]) for line in lines[21 + nb_som + 2 * nb_cel + nb_fac:21 + nb_som + 2 * nb_cel + 2 * nb_fac]]
    fac = [[p - 1 for p in c] for c in tmp_fac]


    mesh = ml.MEDCouplingUMesh("mesh", 3)
    mesh.allocateCells(len(cel))

    for e in range(len(cel)):
        con = []
        for face in cel[e]:
            con.extend(fac[face])
            con.append(-1)
        mesh.insertNextCell(ml.NORM_POLYHED, con[:-1])
    mesh.finishInsertingCells()

    pts = []
    for p in som: pts.extend(p)
    co = ml.DataArrayDouble(pts, len(pts) / 3, 3)
    mesh.setCoords(co)

    mf, d, di, r, ri = mesh.buildDescendingConnectivity()
    mf.setName("mesh")
    mm = ml.MEDFileUMesh.New()
    mm.setMeshAtLevel(0, mesh)
    mm.setMeshAtLevel(-1, mf)

    g = []
    nb_vois = ri.deltaShiftIndex()
    for i in range(mf.getNumberOfCells()):
        if nb_vois[i] == 1: g.append(i)
    grp = ml.DataArrayInt.New(g)
    grp.setName("boundary")
    mm.addGroup(-1, grp)

    mm.write(nom_med, 2)

# maillages 2D
for t, m, (l, d) in (("Cartesian", "cart", (1, 7)), ("Locally_Refined", "ref", (1, 7)), ("Quadrangles", "quad", (1, 7)), ("Triangles", "tri", (1, 6))):
    for n in range(l, d + 1):
        print t, n
        read_typ2("Meshes/2D/{}/mesh_{}_{}.typ2".format(t, m, n), "{}_{}.med".format(t, n))

# maillages 3D
for t, m, (l, d) in (("Hexa", "hexa", (1, 5)), ("Prism", "prism", (1, 4)), ("Tetra", "tetra", (0, 6))):
    for n in range(l, d + 1):
        print t, n
        read_typ3("Meshes/3D/{}/mesh_{}_{}.typ3".format(t, m, n), "{}_{}.med".format(t, n))