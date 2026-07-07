# Szintező algoritmus irányításra
def szintezo(fajl):
    #Fájlbeolvasás
    with open(fajl, "r") as f:
        n = int(f.readline().strip())
        g = list(map(int, f.readline().strip().split(" ")))
        e = [list(map(int, sor.strip().split())) for sor in f]
    eredeti = e.copy()
    
    
    #Csúcsok éllistája bejövő élek szerint, ellista[csúcs] = [szint, [bejövő éleken szomszédok]]
    ellista = {i: [0, []] for i in range(n)}
    for el in e:
        ellista[el[1]][-1].append(el[0])
        
    # Aktív csúcsok meghatározása szintenként aktiv_csucsok[szint] = [aktiv csúcsok a szinten] 
    aktiv_csucsok = {i:[] for i in range(n+1)}
    for csucs in ellista:
        if len(ellista[csucs][1]) > g[csucs]:
            aktiv_csucsok[0].append(csucs)
            
    
    nincs_megoldas = False
    while any(len(v) != 0 for v in aktiv_csucsok.values()):
        #Megkeressük a legmagasabb szinten lévő egyik aktív csúcsot
        for i in reversed(range(n)):
            if aktiv_csucsok[i]:
                aktiv_csucs = aktiv_csucsok[i][0]
                break
        #Ha nincs több él ami bemegy a csúcsba megemeljük
        if len(ellista[aktiv_csucs][1]) == 0:
            ellista[aktiv_csucs][0] += 1
            continue
        
        #Végignézzük, hogy a "beszomszédai" közül van-e olyan, ami egy szinttel lejebb van
        javitas = False
        aktiv_szint = ellista[aktiv_csucs][0]
        for szomszed in ellista[aktiv_csucs][1]:
            if ellista[szomszed][0] ==  aktiv_szint -1:
                #Ha van akkor megfordítjuk az élet
                ellista[aktiv_csucs][1].remove(szomszed)
                ellista[szomszed][1].append(aktiv_csucs)
                
                #Ellenőrizzük, hogy a megfordítás a szomszédot aktívvá teszi-e
                if szomszed not in aktiv_csucsok[aktiv_szint - 1] and len(ellista[szomszed][1]) > g[szomszed]:
                    #Ha igen betesszük a szintjében az aktív csúcsok közé
                    aktiv_csucsok[aktiv_szint - 1].append(szomszed)
                #Ellenőrizzük, hogy az eddigi aktív csúcs pontos lett-e
                if len(ellista[aktiv_csucs][1]) == g[aktiv_csucs]:
                    aktiv_csucsok[aktiv_szint].remove(aktiv_csucs)
                javitas = True
                break

        # Ha nincs javítás, akkor emelünk egyet az aktív csúcs szintjén
        if not javitas:
            ellista[aktiv_csucs][0] += 1
            aktiv_csucsok[aktiv_szint].remove(aktiv_csucs)
            aktiv_csucsok[aktiv_szint+1].append(aktiv_csucs)
        
        #Ha van aktív csúcs az n-dik szinten
        if len(aktiv_csucsok[n]) == 1:
            nincs_megoldas = True
            serto = []
            
            szintek = set()
            for value in ellista.values():
                szintek.add(value[0])
            for i in range(n):
                if i not in szintek:
                    ures = i
                    break
            for csucs, value in ellista.items():
                if value[0] > ures:
                    serto.append(csucs)                   
                
            with open('output.txt', 'w') as f:
                for csucs in serto:
                    f.write(str(csucs)+"\n")
            break
    
    if not nincs_megoldas:
        for el in range(len(eredeti)):
            kezdo = eredeti[el][0]
            veg = eredeti[el][1]
            if kezdo in ellista[veg][1]:
                ellista[veg][1].remove(kezdo)
            else:
                eredeti[el].reverse()
                ellista[kezdo][1].remove(veg)
        with open('output.txt', 'w') as f:
            for el in eredeti:
                f.write(str(el[0])+ ' ' +str(el[1])+'\n')

szintezo('szintezo_proba.txt')