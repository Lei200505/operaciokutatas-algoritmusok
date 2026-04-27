import os
import sympy as sp
from sympy import Matrix, pprint, Rational




def primal_szimplex_lepes(A, c, basis, x_0):
    #Ha kimenet 1, akkor optimális a megoldás 
    #Ha kimenet 0 akkor lépünk egyet, ha kimenet -1 akkor nem korlátos a primál feladat 
    # és visszaad x_0 helyére egy növelőirányban korlátlanul növelő irányt
    kimenet = 1
    
    #y meghatározása
    valtozok_szama = A.shape[1] 
    c_B = sp.Matrix([c[i] for i in basis])
    y = A[:, basis].T.LUsolve(c_B).T
    print(f"y={list(y)}, y*A={list(y*A)} , c={list(c)}")
    
    #ha van olyan oszlop amelyre y*A < c, akkor a primál feladat nem optimális, és lépni kell
    for idx in range(valtozok_szama):
        reduced_cost = c[idx] - (y * A[:, idx])[0]
        if idx not in basis and reduced_cost > 0:
            bekerulo_idx = idx
            print(f"Belépő változó indexe: {bekerulo_idx}\n")
            kimenet = 0
            break
    else:
        return x_0, y, basis, kimenet

    #x' meghatározása
    B = A[:, basis]
    d_B = -B.LUsolve(A[:, bekerulo_idx])
    x_dot = Matrix.zeros(A.shape[1], 1)
    
    for k, i in enumerate(basis):
        x_dot[i, 0] = d_B[k]
    x_dot[bekerulo_idx, 0] = 1
    print(f"x' értéke: {list(x_dot)}")
    
    
    #ha minden x' >= 0, akkor a primál feladat nem korlátos a növelőirányban,
    #és visszaadjuk x helyére a növelőirányban korlátlanul növelő irányt
    if all(val >= 0 for val in x_dot):
        kimenet = -1
        return x_dot, y, basis, kimenet
    
    #különben meghatározzuk a kilépő változót és a lambda értékét, amivel lépünk
    else:
        kimenet = 0
        candidates = [i for i in range(valtozok_szama) if x_dot[i] < 0]
        lamda = min(x_0[i, 0] / -x_dot[i, 0] for i in candidates)
        kikerulo_idx = min(candidates, key=lambda i: x_0[i, 0] / -x_dot[i, 0])
        print(f"Kilépő változó indexe: {kikerulo_idx}, lambda: {lamda}")
    x_0 += lamda * x_dot
    basis.remove(kikerulo_idx)
    basis.append(bekerulo_idx)
    basis.sort()
    
    print(f"Új x értéke: {list(x_0)}, új bázis: {basis}\n")
    
    return x_0, y, basis, kimenet

    

def ketf_szimplex(fajl):
    #Fájlbeolvasás
    with open(fajl) as f:
        n = int(f.readline().strip())
        m = int(f.readline().strip())
        c = Matrix([list(map(Rational, f.readline().split()))]).T
        Ab = Matrix([list(map(Rational, sor.strip().split())) for sor in f]) 
    #Adatellenőrzés
    if Ab.shape != (m,n+1) or c.shape != (n, 1):
        raise ValueError("A megadott adatok mérete nem egyeznek meg az inputtal")

    #Szükséges b >= 0 feltétel biztosítása és mátrix átalakítása
    for i in range(m):
        if Ab[i, n] < 0:
            Ab[i, :] = -Ab[i, :]
    A = Ab[:, :n]
    b = Ab[:, n]
                
    # Összefüggő sorok eltávolítása
    # Egy sort csak akkor törlünk ha Ab-beli része linerárisan függő a már megtartott sorok Ab-beli részeitől
    # (tehát csak a redundáns sorokat töröljük)
    keep = []
    Ab_kept = Matrix.zeros(0, n+1)
    for i in range(Ab.rows):
        Ab_candidate = Ab_kept.col_join(Ab[i, :])
        # csak akkor vesszük fel, ha növeli a rangot
        if Ab_candidate.rank() > Ab_kept.rank():
            Ab_kept = Ab_candidate
            keep.append(i)
    Ab = Ab_kept
    A = Ab[:, :n]
    b = Ab[:, n]
    m = Ab.rows
    
    #------------------------------Első fázis--------------------------------
    A_dot = Ab[:, :n].row_join(Matrix.eye(m))
    c_dot =Matrix.zeros(1,n).row_join(-Matrix.ones(1,m))
    basis = list(range(n, n+m))
    
    #x_0 meghatározása
    x_0 = Matrix.zeros(n+m, 1)
    for i in basis:
        x_0[i] = b[i-n]
    
    
    vege = 0
    print("=" * 40 + "        E L S Ő   F Á Z I S       " + "=" * 40)
    print(f"Bázis: {basis}, x_0: {list(x_0)}\n")
    i = 1
    while vege == 0:
        if max(basis)+1 <= n:
            x_0 = x_0[:n, :]
            basis = [i for i in basis if i <n]
            print(f"Első fázis vége: x_0= {list(x_0)}, bázis: {basis}")
            vege = 1
            break
        elif all(x_0[i] == 0 for i in basis if i >= n):
            basis_0 = [basis[i] for i in range(len(basis)) if basis[i] < n]
            for i in range(len(basis)):
                if basis[i] >= n:
                    for j in range(n):
                        if A_dot[i, j] != 0:
                            basis_0.append(j)
                            break            
            basis = sorted(set(basis_0))    
            
            print(f"Első fázis vége: x_0= {list(x_0)}, bázis: {basis}")
            vege = 1
            break
        print("-" * 40 + f"  Első fázis: {i}. lépés   " + "-" * 40)
        print(f"Aktuális x értéke: {list(x_0)}, bázis: {basis}")
        i += 1
        x_0, y, basis, vege = primal_szimplex_lepes(A_dot, c_dot, basis, x_0)
    #itt ellenőrizzük, hogy az első fázisban talált megoldás megfelel-e a primál feladat feltételeinek,
    # ha nem akkor nincs megoldás
    if vege == 0 or any(x_0[i] > 0 for i in basis if i >= n):
        raise ValueError("Nincs megoldás")
    else:
        x_0 = x_0[:n, :]
    
         
    # (megjegyzés: elvileg az első fázisban a megoldás mindig korlátos lesz, csak az fordulhat elő,
    # hogy a megoldásban maradnak mesterséges változók, amik nem lehetnek pozitívak)
    #------------------------------Második fázis--------------------------------
    print("=" * 40 + "        M Á S O D I K   F Á Z I S       " + "=" * 40)
    vege = 0
    i=1
    while vege == 0:
        print("-" * 40 + f"  Második fázis: {i}. lépés   " + "-" * 40)
        print(f"Aktuális x értéke: {list(x_0)}, bázis: {basis}")
        i +=1
        
        x_0, y, basis, vege = primal_szimplex_lepes(A, c, basis, x_0)
        if vege == -1:
            return f"A primál feladat nem korlátos az {list(x_0)} növelőirányban nem korlátos"
    
    return f"Optimumérték: {int(list(c.T * x_0)[0])}, megfelelő optimális megoldások: x: {list(x_0)}, y: {list(y)}"


file = os.path.join(os.path.dirname(__file__), './tests/input_1.txt')
print(ketf_szimplex(file))