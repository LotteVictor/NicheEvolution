#Copyright (C) 2020 Charlotte Sieger
# 
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

##13.01.2020
##Haploid, annual species in variable environment,  H_t of Normal distribution
#survival of eggs dependent on environment
#(environment variable H_t, individual optimum h_ind, and tolerance g_ind)
#survival of larvae dependent on density
#Trade-off between fertility and niche-width=Th, α= cost of being a generalist
#in a landscape that varies in this environmental parameter


using Distributions #for different distributions and random numbers 
using ArgParse #to parse arguments
using DelimitedFiles #for writedlm - write in .txt and .jl files
using FreqTables #frequency tables -> gives frequencies of entrys in arrays as Dictionary
using StatsBase #statistical analysis

parasource="/home/ubuntu/ParameterDict.jl"
#"/home/charlotte/Simulation_Run/ParameterDict.jl" #gives parameter source -> dictionary!

#constructor for Patches, has arrays for the trait values, its coordinates and the start environmental value as fields
struct TPatch
    v_id::Array{String,1} #array of individual strain ids
    v_hind::Array{Float64, 1} #array of individual niche optimums
    v_gind::Array{Float64, 1} #array of individual niche width = tolerances
    v_tradeoff::Array{Float64, 1} #array of individual trade-off calculated from v_gind
    v_fert::Array{Float64, 1} #array of individual fertilities
    v_a::Array{Float64,1} #array of individual base dispersal probabilities
    v_deltah::Array{Float64,1} #array of individual fitness threshold
    v_dispersed::Array{Bool,1} #did the individual at this position disperse before? 
    row::Int64 #coordinates of the patch: row
    col::Int64 #coordinates of the patch: column
    H_t_start::Float64 #value of the environment from landscapefile in this patch to start climate
end

#function to read in all the necessary arguments, like landscapesource, dispersal type and trade-off strength from a shell skript
function read_arguments()
    s=ArgParseSettings() #setzt s als ArgParseSettings object
    @add_arg_table s begin
        "--landscapesource", "-l" #Names des Argument das eingegeben werden muss in .sh file
            help="Which landscape source file should be used?"
        "--alpha", "-a" #parsing alpha in sh file. alpha is trade-off strength, low alpha -> high trade-off
        help="How strong is the trade-off? Give alpha"
        arg_type=Int64
        default=4
        "--localclimate", "-c" #environment fluctuation locally or globally?
        arg_type=Bool
        default=false
        "--global_dispersal", "-g" #sets dispersal type, toggle between global dispersal (random) and nearest neighbour
        arg_type=Bool
        default=true
        "--habitatemi", "-e" #toggles between habitat depending emigration and stochastic dispersal
        arg_type=Bool
        default=true
        "--minfert" "-f" #gives threshold value for fertility dependent emigration 
        arg_type=Float64
        default=5.0
        "--burninperiod", "-b"  #if true, environment starts to change after a given burn-in period. if set to false, change starts in gen 2
        arg_type=Bool
        default=false
        "--vartrend", "-v" #if true, environmental standard deviation increases with step size specified in parasource
        arg_type=Bool
        default=false
        "--meantrend", "-m" #if true, environmental mean increases with step size specified in parasource
    end
    return parse_args(s)
end

#function to automatically create the filename of the results .txt-file - is also further used to name the files of the frequency table of differents strains and the emmigration and immigration rates
function set_filename(argumentsdict::Dict, par::Dict, n::Int64)
    if argumentsdict["global_dispersal"]==true
        filename_end="_GR"
    elseif (argumentsdict["global_dispersal"]==false) && (argumentsdict["habitatemi"]==false)
        filename_end="_NN"
    elseif argumentsdict["habitatemi"]==true
        filename_end="_HE"
    end
    if argumentsdict["vartrend"]==true
        filename_end=string(filename_end, "_v", string(par["H_t_sd_step"]))
    end
    if argumentsdict["meantrend"]==true
        filename_end=string(filename_end, "_m", string(par["H_t_mean_step"]))
    end
    if argumentsdict["localclimate"]==true
        filename_end=string(filename_end, "_l")
    end
    landscapesource=argumentsdict["landscapesource"]
    filename=string(landscapesource[1:end-18], landscapesource[end-10],landscapesource[end-6], "_",landscapesource[end-4], filename_end, n, "_al2.txt")
    return landscapesource, filename
end

#function to read in the data from the landscape file -> saves it in H_t_start, because the landscape file gives the environmental values to start the simulations with
function read_world(landscapesource::String)
    landscapefile=readdlm(landscapesource)
    numberrows=length(landscapefile[1:end,1])
 #   println(numberrows)
    numbercols=length(landscapefile[1,1:end])
  #  println(numbercols)
    H_t_start=Array{Float64,2}(undef, numberrows,numbercols)
    for i in 1:numberrows
        for j in 1:numbercols
            H_t_start[i,j]=landscapefile[i,j]
        end
    end
    return H_t_start
end

### initialization
#Create the world and initialize the adults 
#Function to initialize Adults as Arrays of length Adultstart
#one individual has several traits:
#hind (individual niche optimum) -- evolves
#gind (individual tolerance /niche width) -- evolves
#tradeoff calculated from gind and the cost to be a generalist: according to Chaianunporn & Hovestadt 2015: Th= exp(-gind²/tcost²)
#a (minimum dispersal probability) -- evolves
#deltah (individual threshold fertility below which individuals disperse) -- evolves
#id (individual id to mark the different strains)
#each trait is one array per patch, individuals are represented by indices
function InitWorld(landscapesource::String, par::Dict, alpha::Int64)
    #println("sigma= ", par.σgind, "mu= ",par.μgind )
    H_t_start=read_world(landscapesource)
    numberrows=length(H_t_start[1:end,1])
    numbercols=length(H_t_start[1,1:end])
    v_hind=Array{Float64,1}(undef, par["Adultstart"])
    v_gind=Array{Float64,1}(undef, par["Adultstart"])
    v_tradeoff=Array{Float64,1}(undef, par["Adultstart"])
    landscape=Array{TPatch,2}(undef, numberrows,numbercols)
    maladapt=Array{Float64,1}(undef, par["Adultstart"])
    v_fert=Array{Float64,1}(undef, par["Adultstart"])
    v_a=Array{Float64,1}(undef, par["Adultstart"])
    v_deltah=Array{Float64,1}(undef, par["Adultstart"])
    v_dispersed=Array{Float64,1}(undef, par["Adultstart"])
    v_id=Array{String,1}(undef, par["Adultstart"])
    for p in 1:numberrows
        for q in 1:numbercols
            v_id=string.(p, "_", q, "_", collect(1:par["Adultstart"]))
            #id is combination of coordinates of the natal patch and the number of the individual
            v_hind=rand(Normal(par["hind_mean_default"],par["hind_std_default"]), par["Adultstart"])
            v_gind=rand(LogNormal(par["μgind_default"],par["σgind_default"]), par["Adultstart"])
            v_tradeoff= exp. (-(0.5 .*v_gind .^2 .*(1/alpha)^2))
            maladapt=((v_hind .- H_t_start[p,q]).^2)./(v_gind.^2) #how well are the individuals adapted to the current environment in their patch
            v_fert =par["R0"] .*v_tradeoff .*exp. (-maladapt)
            v_a=rand(Uniform(0,0.5), par["Adultstart"])
            v_deltah=rand(Uniform(0, 10), par["Adultstart"])
            v_dispersed=fill(false, par["Adultstart"])
            patch= TPatch(v_id, v_hind, v_gind, v_tradeoff, v_fert, v_a, v_deltah,v_dispersed, p, q, H_t_start[p,q])
            landscape[p,q]=patch
        end    
    end
    return landscape, H_t_start
end

#function to create a two-dimensional  array (representing the landscape) where each element contains an array of Float64, with the environmental means for each year.
#-> H_t_landscape is a three-dimensional array
function set_environment(H_t_start::Array{Float64,2}, par::Dict, argumentsdict::Dict)
    numberrows=length(H_t_start[1:end,1])
    numbercols=length(H_t_start[1,1:end])
    H_t_landscape= Array{Array{Float64,1},2}(undef, numberrows, numbercols)
    if argumentsdict["burninperiod"]==true Tini= Int(par["Generations"]/10)
    else Tini = 0
    end
    if argumentsdict["vartrend"]==true
        v_sigma= [fill(par["H_t_sd_default"], Tini); par["H_t_sd_default"] .+ collect((Tini+1):par["Generations"]) .* par["H_t_sd_step"]]
    else
        v_sigma=fill(par["H_t_sd_default"], par["Generations"])    
    end
    if argumentsdict["localclimate"]==true
        for p in 1:numberrows
            for q in 1:numbercols
                if argumentsdict["meantrend"]
                    v_mean= [fill(H_t_start[p,q], Tini); H_t_start[p,q] .+ collect((Tini+1):par["Generations"]) .* par["H_t_mean_step"]]
                else
                    v_mean= fill(H_t_start[p,q], par["Generations"])
                end            
                ran_H= v_mean .+ [rand(Normal(0, v_sigma)) for v_sigma in v_sigma]
                H_t_landscape[p,q] = ran_H
            end
        end
    else
        annualchange=[rand(Normal(0, v_sigma)) for v_sigma in v_sigma]
        for p in 1:numberrows
            for q in 1:numbercols
                if argumentsdict["meantrend"]==true
                    v_mean= [fill(H_t_start[p,q], Tini); H_t_start[p,q] .+ collect((Tini+1):par["Generations"]) .* par["H_t_mean_step"]]
                else
                    v_mean= fill(H_t_start[p,q], par["Generations"])
                end
                ran_H= v_mean .+ annualchange
                H_t_landscape[p,q] = ran_H
            end
        end
    end
    return H_t_landscape
end

#function to get index of target patch under global dispersal 
function global_patch(numbercols::Int64,numberrows::Int64,p::Int64,q::Int64) 
    targetrow = rand(1:numberrows)
    targetcol = rand(1:numbercols)
    if targetrow == p && targetcol==q
        targetcol = rand(1:numbercols)
        targetrow = rand(1:numberrows)        
    end
    if targetrow > numberrows targetrow = 1 end
    if targetcol > numbercols targetcol = 1  end
    if targetrow < 1 targetrow = numberrows  end
    if targetcol < 1 targetcol = numbercols end
    return targetrow, targetcol
end

#function for nearest neighbour dispersal, p and q are row and column in landscape
function nearest_patch(numberrows::Int64, numbercols::Int64, p::Int64, q::Int64)
    targetrow = p+ rand(-1:1)
    targetcol=q+ rand(-1:1)
    if targetrow == p && targetcol==q
        targetrow = p+ rand(-1:1)
        targetcol = q+ rand(-1:1)
    end
    if targetrow > numberrows targetrow = 1 end
    if targetcol > numbercols targetcol = 1 end
    if targetrow < 1 targetrow = numberrows end
    if targetcol < 1 targetcol = numbercols end
    return targetrow, targetcol
end

#function to CHANGE the landscape! Changes the position of the individuals in the landscape, or kills them (death through dispersal), also counts how many individuals emigrated and immigrated for each patch and returns those results in two two-dimensional arrays 
function dispersal!(landscape::Array{TPatch,2}, argumentsdict::Dict, par::Dict)
    numberrows=length(landscape[1:end,1])
    numbercols=length(landscape[1,1:end])
    dispprob=Array{Array{Float64,1},2}(undef, numberrows, numbercols)
    n_disp=zeros(Int64, numberrows, numbercols)  #number of emigrated individuals
    n_imm=zeros(Int64, numberrows, numbercols) #number of immigrated individuals
    for p in 1:numberrows
        for q in 1:numbercols
            absdiff= abs. (landscape[p,q].H_t_start .- landscape[p,q].v_hind)
            probvec=Array{Float64, 1}(undef, length(landscape[p,q].v_hind))
            if (argumentsdict["habitatemi"]==false) #if habitat independent dependent dispersal
                for j in 1: length(landscape[p,q].v_hind)
                    probvec[j]=landscape[p,q].v_a[j] #This is necessary to not have pointer references 
                end #dispersal probability is the individual's base probability
            elseif (argumentsdict["habitatemi"]==true) #if habitat dependent dispersal
                for j in 1: length(landscape[p,q].v_hind)
                    if (landscape[p,q].v_deltah[j] < absdiff[j]) #and if difference between hind and environment is below the threshold
                        probvec[j]=1 #dispersal probability is one
                    elseif landscape[p,q].v_deltah[j]> absdiff[j] #but if difference between hind and environment is above the threshold
                        probvec[j]=landscape[p,q].v_a[j] #dispersal probability is base probability
                    end
                end
            end
            dispprob[p,q]=probvec
        end
    end
    
    for p in 1:numberrows
        for q in 1:numbercols
            i=1
            while i<length(landscape[p,q].v_hind) #looping through all individuals in that patch
                if (argumentsdict["habitatemi"]==false) #if habitat independent dispersal
                    if argumentsdict["global_dispersal"]==true #and global dispersal
                        newrow, newcol= global_patch(numberrows,numbercols, p, q) #the new patch is global
                    elseif argumentsdict["global_dispersal"]==false #and nearest neighbour dispersal
                        newrow, newcol= nearest_patch(numberrows, numbercols, p, q) #the new patch is in the neighboorhood
                    end
                elseif (argumentsdict["habitatemi"]==true) #if habitat dependent dispersal
                    if argumentsdict["global_dispersal"]==true  #and global dispersal
                    newrow, newcol= global_patch(numberrows,numbercols, p, q) #the new patch is global
                    elseif argumentsdict["global_dispersal"]==false #and nearest neighbour dispersal
                    newrow, newcol= nearest_patch(numberrows, numbercols, p, q) #the new patch is in the neighboorhood
                    end
                end
                if (landscape[p,q].v_dispersed[i] == false) &&(dispprob[p,q][i]<rand()) &&(rand()>par["dispmort"]) #if  individual survives dispersal
                    n_disp[p,q]+=1 #number of emigrated individuals increases
                    n_imm[newrow,newcol]+=1 #number of immigrated individuals increases
                    #put dispersed individual into new patch
                    push!(landscape[newrow,newcol].v_id, landscape[p,q].v_id[i]) 
                    push!(landscape[newrow,newcol].v_hind, landscape[p,q].v_hind[i]) 
                    push!(landscape[newrow,newcol].v_gind, landscape[p,q].v_gind[i])
                    push!(landscape[newrow,newcol].v_tradeoff, landscape[p,q].v_tradeoff[i])
                    push!(landscape[newrow,newcol].v_a, landscape[p,q].v_a[i])
                    push!(landscape[newrow,newcol].v_deltah, landscape[p,q].v_deltah[i])
                    push!(landscape[newrow, newcol].v_dispersed, true)
                    push!(dispprob[newrow,newcol], dispprob[p,q][i])
                    #remove individual from natal patch
                    splice!(landscape[p,q].v_id, i)
                    splice!(landscape[p,q].v_hind, i)
                    splice!(landscape[p,q].v_gind, i)
                    splice!(landscape[p,q].v_tradeoff, i)
                    splice!(landscape[p,q].v_a, i)
                    splice!(landscape[p,q].v_deltah, i)
                    splice!(landscape[p,q].v_dispersed, i)
                    splice!(dispprob[p,q], i) #also remove dispersal probability from array, otherwise the indices will get messed up!
                    i -=1
                elseif (landscape[p,q].v_dispersed[i] == false) &&(rand()<dispprob[p,q][i]) &&(rand()<par["dispmort"])
                    #if dispersal but death during migration
                    #remove individual from natal patch
                    n_disp[p,q]+=1 #number of emigrated individuals increases, but not of immigrated individuals -> they die before that!
                    splice!(landscape[p,q].v_id, i)
                    splice!(landscape[p,q].v_hind, i) 
                    splice!(landscape[p,q].v_gind, i)
                    splice!(landscape[p,q].v_tradeoff, i)
                    splice!(landscape[p,q].v_a, i)
                    splice!(landscape[p,q].v_deltah, i)
                    splice!(landscape[p,q].v_dispersed, i)
                    splice!(dispprob[p,q], i) #also remove dispersal probability from array, otherwise the indices will get messed up!
                   # println(length(dispprob[p,q]))
                    i -=1
                end
                i += 1
            end
        end
    end
return n_disp, n_imm
end

#function for each generation, including dispersal, reproduction and mutation
function NextGenFu(landscape::Array{TPatch,2}, H_t_landscape::Array{Array{Float64,1},2}, currGen::Int64, argumentsdict::Dict, par::Dict)
    numberrows=length(landscape[1:end,1])
    numbercols=length(landscape[1,1:end])
    for p in 1:numberrows
        for q in 1:numbercols
            if length(landscape[p,q].v_hind)>0
                #initializing all needed arrays
                Offspring=Array{Int64,1}(undef, length(landscape[p,q].v_hind))
                maladapt=Array{Float64,1}(undef, length(landscape[p,q].v_hind))
                fert=Array{Float64,1}(undef, length(landscape[p,q].v_hind))
                maladapt=((landscape[p,q].v_hind .- H_t_landscape[p,q][currGen]).^2)./(landscape[p,q].v_gind.^2) #how well are the individuals adapted to the current environment in their patch?
                fert .=par["R0"] .*landscape[p,q].v_tradeoff .*exp. (-maladapt) #how fertile are they?
                Offspring=[rand(Poisson(fert)) for fert in fert] # random number of eggs for each indidividual with the fertility of each offspring as the mean of the Poisson distribution
                nOff=sum(Offspring) #total number of offspring in the whole population
                psurv=1/(1+par["a"]*nOff) #survival probability
                survOff=[rand(Binomial(Offspring, psurv)) for Offspring in Offspring] #draw the actual number of offspring per individual
                #empty arrays to push new values in
                v_id2=Array{String,1}() 
                v_hind2=Array{Float64,1}() 
                v_gind2=Array{Float64,1}()
                v_fert2=Array{Float64,1}()
                v_a2=Array{Float64,1}()
                v_deltah2=Array{Float64,1}()
                v_dispersed=fill(false, sum(survOff))
                #cycle through all individuals and repeat them as often as they had offspring
                for i in 1:length(survOff)
                    append!(v_id2, fill(landscape[p,q].v_id[i], survOff[i]))
                    append!(v_hind2, fill(landscape[p,q].v_hind[i], survOff[i]))
                    append!(v_gind2, fill(landscape[p,q].v_gind[i], survOff[i]))
                    append!(v_fert2, fill(fert[i], survOff[i]))
                    append!(v_a2, fill(landscape[p,q].v_a[i], survOff[i]))
                    append!(v_deltah2, fill(landscape[p,q].v_deltah[i], survOff[i]))
                end
                # mutations
                v_hind=v_hind2 .+rand(Normal(0, par["mut_h"]), length(v_hind2)) # mutate optimum trait
                v_gind=v_gind2 .*rand(Uniform(1-par["mut_g"], 1+par["mut_g"]), length(v_gind2)) # mutate tolerance trait
                v_tradeoff=exp. (-(0.5 .*v_gind .^2 .*(1/argumentsdict["alpha"])^2))
                v_a=v_a2 .+rand(Normal(0, 0.001), length(v_a2)) #mutate the base dispersal probability
                v_deltah=v_deltah2 .+rand(Normal(0, 0.001), length(v_deltah2)) #mutate the fertility threshold
                patch= TPatch(v_id2, v_hind, v_gind, v_tradeoff, v_fert2, v_a, v_deltah, v_dispersed, p, q, H_t_landscape[p,q][currGen]) #create a new patch with those arrays (changing the existing patch does not work unless it's a mutable struct, which requires lots of calculation time -> making anew patch and puttin git in the landscape is faster!
                landscape[p,q]=patch #put the new patch on the old position -> replacing the old patch
            end
        end
    end
    return landscape
end

#function to store the data in a text file with a given frequency
function analyze(landscape::Array{TPatch,2}, H_t_start::Array{Float64,2}, H_t_landscape::Array{Array{Float64, 1},2}, n_disp::Array{Int64,2}, n_imm::Array{Int64,2}, filename::String, currGen::Int64)
    numberrows=length(landscape[1:end,1])
    numbercols=length(landscape[1,1:end])
    for p in 1:numberrows
        for q in 1:numbercols
            freqtable_id= Dict(collect(sort(countmap(landscape[p,q].v_id)))[(end-9):end])
            open(filename, "a") do IO
            writedlm(IO, [currGen landscape[p,q].row landscape[p,q].col H_t_landscape[p,q][currGen] H_t_start[p,q] length(landscape[p,q].v_hind) mean(landscape[p,q].v_hind) std(landscape[p,q].v_hind) mean(landscape[p,q].v_gind) std(landscape[p,q].v_gind) mean(landscape[p,q].v_fert) mean(landscape[p,q].v_a) mean(landscape[p,q].v_deltah)])
            end #write to general results file
            open(string(filename[1:end-4], "freqtable.jl"), "a") do IO
            writedlm(IO, [currGen landscape[p,q].row landscape[p,q].col freqtable_id])
            end #write to frequency of strains results file
            open(string(filename[1:end-4], "dispersalrates.txt"), "a") do IO
            writedlm(IO, [currGen landscape[p,q].row landscape[p,q].col length(landscape[p,q].v_hind) n_disp[p,q] n_imm[p,q]])
            end #write to dispersal rates results file
        end
    end
end

#One funciton to rule them all! This finction actually runs the simulation by calling the other functions in the correct order
function Simulation_Run(parasource::String, n::Int64)
    println("Now starting Simulation") #So we now, the simulation runs
    par=include(parasource)
    argumentsdict=read_arguments()
    colnames=["generation" "row" "col" "curr_env" "mean_env" "N" "meanh" "stdh" "meang" "stdg" "meanfert" "mean_dispprob" "mean_minfert"] #the column names for the results file
    landscapesource,filename =set_filename(argumentsdict, par, n) #setting the filename and reading the landscapesource
    landscape, H_t_start =InitWorld(landscapesource, par, argumentsdict["alpha"])#initializing the world
    H_t_landscape=set_environment(H_t_start, par, argumentsdict)#creating climate
    println("Environment set")
    writedlm(filename, colnames) #create the results file
    writedlm(string(filename[1:end-4], "freqtable.jl"), ["Generation" "Row" "Col" "Freq_ID"]) #create the frequency table file -> this is a -jl file, because the freq table is a dictionary!
    writedlm(string(filename[1:end-4], "dispersalrates.txt"), ["Generation" "Row" "Col" "N" "Emmigrated" "Immigrated"]) #create the dispersal rates file
    for i in 2:par["Generations"]
        n_disp, n_imm = dispersal!(landscape, argumentsdict, par) #dispersal before reproduction
        if i%10==0 #every tenth generation
            analyze(landscape,H_t_start, H_t_landscape, n_disp, n_imm, filename, i) #append all three results file with the current values
           # println(i, " Generations complete")
        end
        landscape =NextGenFu(landscape, H_t_landscape, i, argumentsdict, par) #next generation
    end
end

@time for j in 1: 3
    Simulation_Run(parasource, j) #actually run the simulation and time it!
end
