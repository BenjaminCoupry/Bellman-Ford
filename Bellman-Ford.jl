l = Dict(("A","B")=>3,("A","E")=>5,("E","B")=>1,
("B","C")=>4,("C","D")=>2,("E","D")=>9,("D","F")=>3)
l1 = Dict(("A","B")=>0.3,("A","E")=>0.5,("E","B")=>0.1,
("B","C")=>0.4,("C","D")=>-0.2,("E","D")=>0.9,("D","F")=>0.3,("D","B")=>-0.3)
noms = ["A","B","C","D","E","F"]

disjonctions = [[("D","C")=>10,("C","D")=>20],[("F","C")=>11,("F","D")=>21]]
choix = [2,2]

taches = Dict("A"=>(2,[],[]),"B"=>(3,["A"],[]),"C"=>(1,["A"],[]),"D"=>(4,["B","C"],[]),"E"=>(1,["C"],[]))
travaux = [[(6,["M1"]),(7,["M3"])],[(3,["M1"]),(5,["M2"]),(1,["M3"])]]

function getPoids(i,j, Liaisons, Noms)
    return Liaisons[(Noms[i],Noms[j])]
end
function arcExiste(i,j, Liaisons, Noms)
    return haskey(Liaisons,(Noms[i],Noms[j]))
end
function getIndice(nom, Liaisons, Noms)
    nb = length(Noms)
    for i in 1:nb
        if Noms[i] == nom
            return i
        end
    end
end

function getSuccesseursTache(nom_tache,Taches)
    res = []
    for ntache in keys(Taches)
        tache = Taches[ntache]
        pred = tache[2]
        if nom_tache in pred
            push!(res,ntache)
        end
    end
    return res
end
#Taches est un dict de taches
#une tache est nom=>(duree,[] des noms des predecesseurs,[] des noms des ressources utilisees)
function Taches2Transitions(Taches)
    noms = ["Debut","Fin"]
    liaisons = Dict()
    usageRessources = Dict()
    for ntache in keys(Taches)
        push!(noms,ntache)
        tache = Taches[ntache]
        succ = getSuccesseursTache(ntache,Taches)
        pred = tache[2]
        duree = tache[1]
        ress = tache[3]
        if length(succ) ==0
            liaisons[(ntache,"Fin")] = duree
        end
        if length(pred) ==0
            liaisons[("Debut",ntache)] = 0
        end
        for successeur in succ
            liaisons[(ntache,successeur)] = duree
        end
        for ressource in ress
            if ressource != ""
                if haskey(usageRessources,ressource)
                    push!(usageRessources[ressource],(ntache,duree))
                else
                    usageRessources[ressource]=[(ntache,duree)]
                end
            end
        end
    end
    djctions = usageRessource2Disjonctions(usageRessources)
    return noms,liaisons,djctions
end

#usageRessources est un dict, les clefs sont les noms des ressources, les valeurs sont la liste des  (nom,temps) des taches qui utilisent cette ressource. 
function usageRessource2Disjonctions(UsageRessources)
    Disjonctions = []
    for ressource in keys(UsageRessources)
        usage = UsageRessources[ressource]
        n = length(usage)
        for moi in 1:n 
            nom_moi = usage[moi][1]
            temps_moi = usage[moi][2]
            for lui in moi+1:n 
                nom_lui = usage[lui][1]
                temps_lui = usage[lui][2]
                disjonct = [(nom_moi,nom_lui)=>temps_moi,(nom_lui,nom_moi)=>temps_lui]
                push!(Disjonctions,disjonct)
            end
        end
    end
    return Disjonctions
end
# travaux est une liste de travail
#un travail est une liste de taches
#une tache est (temps,[] de nom_ressource utilisees), mettre "" si pas de ressource

function travaux2Transitions(Travaux)
    noms = ["Debut","Fin"]
    liaisons = Dict()
    usageRessources = Dict()
    tr=1
    for travail in Travaux
        ta=1
        precedent = "Debut"
        temps_preced=0
        for tache in travail
            nom = "travail"*string(tr)*"_tache"*string(ta)
            push!(noms,nom)
            temps = tache[1]
            ressources = tache[2]
            for ressource in ressources
                if ressource != ""
                    if haskey(usageRessources,ressource)
                        push!(usageRessources[ressource],(nom,temps))
                    else
                        usageRessources[ressource]=[(nom,temps)]
                    end
                end
            end
            liaisons[(precedent,nom)] = temps_preced
            temps_preced=temps
            precedent = nom
            ta+=1
        end
        liaisons[(precedent,"Fin")] = temps_preced
        tr+=1
    end
    djctions = usageRessource2Disjonctions(usageRessources)
    return noms,liaisons,djctions
end

#ajoute aux liaisons les disjonctions données en prenant pour la ieme disjonction l'option spécifiée par le ieme element de choix,choix[i]=0 pour ne rien choisir
function ajouterLiaisonsDisjonctions(Liaisons, Disjonctions, Noms, Choix)

    liaisons_etendues = deepcopy(Liaisons)
    for i in 1:length(Choix)
        possibilites = Disjonctions[i]
        choix_local = Choix[i]
        if choix_local != 0 
            kv = possibilites[choix_local]
            k,v=kv
            if haskey(liaisons_etendues,k)
                liaisons_etendues[k]=max(v,liaisons_etendues[k])
            else
                liaisons_etendues[k]=v
            end
        end
    end
    return liaisons_etendues,noms
end

function bellmanFord(Start, minmax, modecalc,iterations, Liaisons, Noms)
    nbSommets = length(Noms)
    nbLiaisons = length(keys(Liaisons))
    if iterations == -1
        nbIt = nbSommets
    else 
        nbIt = iterations
    end
    start = getIndice(Start, Liaisons, Noms)
    #Preparation
    Resultats = zeros(nbIt,nbSommets)
    Passages = Array{String,2}(undef, nbIt, nbSommets)

    # pires valeurs possibles
    for sommet in 1:(nbSommets)
        if minmax == "min"
            Resultats[1,sommet] = Inf
        elseif minmax == "max"
            Resultats[1,sommet] = -Inf
        end
        Passages[1,sommet] = "_"
    end

    #Elements neutres pour l'operation
    if modecalc == "somme"
        Resultats[1,start] = 0
    elseif modecalc == "plusFaible"
        Resultats[1,start] = Inf
    elseif modecalc == "plusFort"
        Resultats[1,start] = -Inf
    elseif modecalc == "produit"
        Resultats[1,start] = 1
    end
    
    #Calcul
    for iteration in 2:nbIt
        for sommet in 1:nbSommets
            #Distance par défaut
            dcalc = Resultats[iteration-1,sommet]
            Passages[iteration,sommet] = Passages[iteration-1,sommet]
            for candidat in 1:nbSommets
                #La transition est elle validee
                #distance entre l'origine et le sommet de candidat tésté
                start2cand = Resultats[iteration-1,candidat]
                produitApprouved = ( modecalc != "produit" || (modecalc == "produit" && abs(start2cand) != Inf))
                valide = arcExiste(candidat,sommet, Liaisons, Noms) && produitApprouved
                if valide
                    #distance entre le sommet de candidat tésté et le sommet en cours de calcul
                    poidsArc = getPoids(candidat,sommet, Liaisons, Noms)
                    # Calcul de la fonction de cout du candidat
                    if modecalc =="plusFaible"
                        dcandidat = min(start2cand,poidsArc)
                    elseif modecalc =="plusFort"
                        dcandidat = max(start2cand,poidsArc)
                    elseif modecalc == "somme"
                        dcandidat = start2cand+poidsArc
                    elseif modecalc == "produit"
                        dcandidat = start2cand*poidsArc
                    end
                    #on garde la minimale ou maximale en fonction des cas
                    condmin = (dcandidat<dcalc)
                    condmax = (dcandidat>dcalc)
                    if((condmin && (minmax == "min")) || (condmax &&(minmax == "max")))
                        dcalc = dcandidat
                        Passages[iteration,sommet] = Noms[candidat]
                    end
                end
            end

            Resultats[iteration,sommet] = dcalc
        end
    end
    display(Resultats)
    display(Passages)
    return Resultats,Passages
end

#convertit une table de transition en dictionaire des transitions
function tableTransition2Dict(Table)
    n = size(Table)[1]
    noms = []
    transitions = Dict()
    for i in 1:n
        push!(noms,string(i))
    end
    for i in 1:n
        for j in 1:n
            if Table[i,j] != Missing
                transitions[(string(i),string(j))] = Table[i,j]
            end
        end
    end
    return noms,transitions
end

# depart, arrivee, max/min = chemin maximisant ou minimisant
# plusFaible/plusFort/somme = mode de calcul de l'effet des arretes
# -1 pour nombre d'iterations par defaut, n pour n iterations
#noms la liste des noms des Sommets
#l2 le dict des transitions, clefs (nomDepart,nomArrivee), valeurs poids des arcs
#display(chemin("A","F","max","produit",-1, l2, noms))
function chemin(Start, Finish, minmax,modecalc,nbIterations, Liaisons, Noms)
    try
        dist,noms = bellmanFord(Start, minmax,modecalc,nbIterations, Liaisons, Noms)
        it = size(noms,1)
        start = getIndice(Start, Liaisons, Noms)
        finish = getIndice(Finish, Liaisons, Noms)
        result = [Finish]
        cible=finish
        while cible != start
            cible = getIndice(noms[it,cible], Liaisons, Noms)
            ncible = Noms[cible]
            push!(result,ncible)
            it-=1
        end
        push!(result,"Cout : " * string(dist[size(noms,1),finish]))
        return reverse(result)
    catch err
        return ["Pas de solution"]
    end
end

# depart, arrivee, max/min = chemin maximisant ou minimisant
# plusFaible/plusFort/somme = mode de calcul de l'effet des arretes
# -1 pour nombre d'iterations par defaut, n pour n iterations
#l2,noms=ajouterLiaisonsDisjonctions(l1, disjonctions, noms, choix)
#display(chemin("A","F","max","produit",-1, l2, noms))

#noms,liaisons,djctions = travaux2Transitions(travaux)
noms,liaisons,djctions = Taches2Transitions(taches)
display(noms)
display(liaisons)
display(djctions)
