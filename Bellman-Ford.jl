Liaisons = Dict(("A","B")=>3,("A","E")=>5,("E","B")=>-1,
("B","C")=>4,("C","D")=>2,("E","D")=>9,("D","F")=>3)

Noms = ["A","B","C","D","E","F"]

function getPoids(i,j)
    return Liaisons[(Noms[i],Noms[j])]
end
function arcExiste(i,j)
    return haskey(Liaisons,(Noms[i],Noms[j]))
end
function getIndice(nom)
    nb = length(Noms)
    for i in 1:nb
        if Noms[i] == nom
            return i
        end
    end
end

function bellmanFord(Start)
    nbSommets = length(Noms)
    nbLiaisons = length(keys(Liaisons))
    nbIt = nbSommets-1
    start = getIndice(Start)
    #Preparation
    Resultats = zeros(nbIt,nbSommets)
    Passages = Array{String,2}(undef, nbIt, nbSommets)

    for sommet in 1:(nbSommets)
        Resultats[1,sommet] = Inf
        Passages[1,sommet] = "_"
    end
    Resultats[1,start] = 0
    
    #Calcul
    for iteration in 2:nbIt
        for sommet in 1:nbSommets
            #Distance par défaut
            dcalc = Resultats[iteration-1,sommet]
            Passages[iteration,sommet] = Passages[iteration-1,sommet]
            for candidat in 1:nbSommets
                if arcExiste(candidat,sommet)
                    #distance entre l'origine et le sommet de candidat tésté
                    start2cand = Resultats[iteration-1,candidat]
                    #distance entre le sommet de candidat tésté et le sommet en cours de calcul
                    poidsArc = getPoids(candidat,sommet)
                    #on garde la minimale
                    dcandidat = start2cand+poidsArc
                    if(dcandidat<dcalc)
                        dcalc = dcandidat
                        Passages[iteration,sommet] = Noms[candidat]
                    end
                end
            end
            Resultats[iteration,sommet] = dcalc
        end
    end

    return Resultats,Passages
end

function chemin(Start, Finish)
    dist,noms = bellmanFord(Start)
    it = size(noms,1)
    start = getIndice(Start)
    finish = getIndice(Finish)
    result = [Finish]

    cible=finish
    while cible != start
        cible = getIndice(noms[it,cible])
        ncible = Noms[cible]
        push!(result,ncible)
        it-=1
    end
    return reverse(result)
end


chemin("A","C")