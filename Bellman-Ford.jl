l = Dict(("A","B")=>3,("A","E")=>5,("E","B")=>1,
("B","C")=>4,("C","D")=>2,("E","D")=>9,("D","F")=>3)
l1 = Dict(("A","B")=>0.3,("A","E")=>0.5,("E","B")=>0.1,
("B","C")=>0.4,("C","D")=>-0.2,("E","D")=>0.9,("D","F")=>0.3,("D","B")=>-0.3)
noms = ["A","B","C","D","E","F"]


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
display(chemin("A","F","max","produit",-1, l1, noms))