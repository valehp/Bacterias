# -*- coding: utf-8 -*-
from Bacterias_v1 import *
import pandas as pd
import openpyxl as oex
import argparse
import os
import time
import numpy as np
import Lexico
import pickle


def data_to_excel(DATOS, df_bests, excel_name='BE_lexico_%.xlsx', num_prueba=1, curso="", df_time=""):
    pandita = pd.DataFrame()
    P, I, C, CL, T, TL, M, ML, LM, LAST_F, GLOBAL_F, IT_GLOBAL, TIME, pob_actual, KP, MULT, EC = [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []

    # (IT, POBLACION, PC, PCL, PT, PTL, PM, PML, PL, min(mi_bacteria.F), mi_bacteria.best_fitness, fitness, it, len(mi_bacteria.P), end-start, estanc_cont, KP, multiplicador) 
    #   0     1       2   3    4    5    6   7   8           9                     10                 11    12          13            14             15     16      17                                  
    for i in range(len(DATOS)):
        P.append(DATOS[i][1])
        I.append(DATOS[i][0])
        C.append(DATOS[i][2])
        CL.append(DATOS[i][3])
        T.append(DATOS[i][4])
        TL.append(DATOS[i][5])
        M.append(DATOS[i][6])
        ML.append(DATOS[i][7])
        LM.append(DATOS[i][8])
        LAST_F.append(DATOS[i][9])
        GLOBAL_F.append(DATOS[i][10])
        IT_GLOBAL.append(DATOS[i][12])
        #pob_actual.append(DATOS[i][13])
        TIME.append(DATOS[i][14])
        KP.append(DATOS[i][16])
        #MULT.append(DATOS[i][17])
        EC.append(DATOS[i][15])
        
    pandita["Población"] = P
    pandita["Iteraciones"] = I
    pandita["C"] = C
    pandita["CL"] = CL
    pandita["T"] = T
    pandita["TL"] = TL
    pandita["M"] = M
    pandita["ML"] = ML
    pandita["Liberar material"] = LM
    pandita["iters para eliminar"] = KP
    #pandita["Multiplicador"] = MULT
    pandita["Mejor F última it"] = LAST_F
    pandita["Mejor F global"] = GLOBAL_F
    pandita["Iteración mejor f global"] = IT_GLOBAL
    #pandita["Población actual"] = pob_actual
    pandita["Tiempo"] = TIME

    d = excel_name.split('/')
    excel_name = d[-1]
    i = 0
    while i < len(d):
        if d[i] == '.' or '.' in d[i]:
            d.pop(i)
            continue
        i += 1

    d = '/'.join(d)
    path = os.getcwd()
    if not d in os.getcwd(): os.chdir(d)

    files = [ file for file in os.listdir(os.getcwd()) ]
    write_mode = 'a' if(excel_name in files) else 'w'
    sheetname_time = curso + "(time)"

    if curso:
        sheet_name = curso
        sheet_name_bests = curso + "(f_it)"
    else:
        sheet_name = "Prueba{} ({})".format(num_prueba, P[0] if(args.change_pob) else I[0])
        sheet_name_bests = "FxIT{} ({})".format(args.num_prueba, P[0] if(args.change_pob) else I[0])



    writer = pd.ExcelWriter(excel_name, engine='openpyxl', mode=write_mode)
    pandita.to_excel(writer, sheet_name=sheet_name )
    df_bests.to_excel(writer, sheet_name=sheet_name_bests)
    df_time.to_excel(writer, sheet_name=sheetname_time)
    writer.save()
    os.chdir(path)


def crear_grupos(g, max_len):
    div = len(g) // max_len
    mod = len(g) % max_len
    grupos = []
    indx = []
    if mod != 0:
        min_len = max_len - 1
        dif = min_len - mod
        grupo_menor = 1 + dif
        grupo_mayor = div - dif
        lens  = [ max_len for i in range(grupo_mayor) ]
        lens += [ min_len for i in range(grupo_menor) ]
        cont, j = 0, 0 # contador de grupos
        for i in range(len(g)):
            if i != 0 and (i-j) % lens[cont] == 0:
                grupos.append( g[j:i] )
                indx.append(j)
                j = i
                cont += 1
            elif i == len(g)-1:
                grupos.append( g[j:] )
                indx.append(j)
    else: 
        j = 0
        for i in range(len(g)):
            if i != 0 and i % max_len == 0:
                grupos.append( g[j:i] )
                indx.append(j)
                j = i
            elif i == len(g)-1:
                grupos.append( g[j:] )
                indx.append(j)
    return grupos, indx

def ret_groups(g, indices):
    all_groups = []
    for index in range(1, len(indices)+1):
        if index >= len(indices):
            i = indices[index-1]
            j = len(g)
        else: 
            i = indices[index-1]
            j = indices[index]
        group = g[i:j]
        all_groups.append(group)
    return all_groups


def UsarBacteria(EVALUACIONES, NOTAS, GRUPOS, ALUMNOS, REPETICIONES, IT, args, save=True):
    """
    - EVALUACIONES:
    - NOTAS:
    - GRUPOS:
    - ALUMNOS:
    - REPETICIONES:
    - IT:
    - args: 
    - save
    """
    DATOS = []
    mejor_fitness = lambda f1, f2: (f1,False) if(f1<=f2) else (f2,True)
    function = None
    if args.fitness == 'lexico': function = Lexico.F
    mejores = []
    if args.elitista: mi_bacteria = BacteriaElitista(EVALUACIONES, NOTAS, GRUPOS, ALUMNOS, args.pob_size, args.PC, args.PCL, args.PT, args.PTL, args.PM, args.PML, args.PL, fitness=function, max_group=args.max_group)
    else: mi_bacteria = Bacteria(EVALUACIONES, NOTAS, GRUPOS, ALUMNOS, args.pob_size, args.PC, args.PCL, args.PT, args.PTL, args.PM, args.PML, args.PL, fitness=function, max_group=args.max_group)

    sol_x_its = []
    df_fits_x_its = pd.DataFrame()
    df_time = pd.DataFrame()
# ----- Empezamos a iterar ----- #
    print( "-"*60 )
    for i in range(REPETICIONES):
        print("rep ", i)
        aux_time = []
        mi_bacteria.reset()
        start = time.time()
        mi_bacteria.generar_poblacion()
        fitness = float('inf')
        it = -1
        fit_cont, estanc_cont, last_kill = 0, 0, 0
        tmp = []
        sol_x_its = []

        for j in range(IT):
            
            if args.elitista:
                if not args.stuck and abs(last_kill - j) == args.iterations_to_kill:
                    last_kill = j
                    mi_bacteria.kill(distint = (True if args.kill == 'f' else False) )
                elif args.stuck and estanc_cont == args.iterations_to_kill:
                    estanc_cont = 0
                    mi_bacteria.kill(distint = (True if args.kill == 'f' else False) )
                   
            mi_bacteria.global_fitness() 
            mi_bacteria.crear_antib()    
            mi_bacteria.clasification()
            tmp = mejor_fitness(fitness, min(mi_bacteria.F))
            if tmp[1]:
                fitness = tmp[0]   
                it = j
                fit_cont = 0

            #print("VGL", end=" -> ")    
            # Variación genética
            mi_bacteria.Conjugacion()
            mi_bacteria.Transformacion()
            mi_bacteria.Mutacion()

            # Reclasificación
            mi_bacteria.global_fitness()
            mi_bacteria.clasification()
            tmp = mejor_fitness(fitness, min(mi_bacteria.F))
            if tmp[1]:
                fitness = tmp[0]
                it = j
                fit_cont = 0

            # Aplicar antibiótico
            mi_bacteria.Aplicar_antibiotico()
            #print("regenerando", end=" -> ")    
            mi_bacteria.Regen()

            # Revisar fitness para ver si se llegó a un estancamiento o no
            a = mi_bacteria.check_fitness(j)
            p_aux = len(mi_bacteria.P)
            if a and args.stuck: estanc_cont += 1
            fit_cont += 1
            sol_x_its.append( min(mi_bacteria.F) )
            #print("LISTO!")    
            aux_time.append( time.time()-start )
            
        end = time.time()
        df_fits_x_its["rep_"+str(i)] = sol_x_its
        df_time["rep_"+str(i)] = aux_time
        tmp_apnd =  (IT, args.pob_size, args.PC, args.PCL, args.PT, args.PTL, args.PM, args.PML, args.PL, min(mi_bacteria.F), mi_bacteria.best_fitness, fitness, it, len(mi_bacteria.P), end-start, estanc_cont, args.iterations_to_kill, 0) 
        DATOS.append(tmp_apnd)
        #print( IT, args.pob_size, args.PC, args.PCL, args.PT, args.PTL, args.PM, args.PML, args.PL, mi_bacteria.best_fitness, it, len(mi_bacteria.P), end-start, estanc_cont, args.iterations_to_kill, 0 )
        mejores.append( mi_bacteria.get_BestSolution() )
    mejores.sort(key=lambda val: val[1] )
    print( vars(args), "\n\t LISTOO!!\n" )
    if save: data_to_excel(DATOS, df_fits_x_its, args.save_fit, num_prueba=args.num_prueba, curso=args.path.split('.')[1].split('/')[-1], df_time=df_time)
    return mejores[0]


def Guardar(excel_name, m, F, DATOS, iteraciones):
    """
    - excel_name: Nombre del excel donde se guardarán los datos
    - write_mode: Modo escritura -> w: no existe el archivo ; a: el archivo existe
    - m: Mejor Bacteria -> Representación de todos los grupos
    - DATOS: Diccionario con los datos de las frecuencias de cada alumno (sacado del .pkl)
    - iteraciones: Número de iteraciones con las que se utilizó la bacteria
    """
    d = excel_name.split('/')
    excel_name = d[-1]
    i = 0
    while i < len(d):
        if d[i] == '.' or '.' in d[i]:
            d.pop(i)
            continue
        i += 1

    d = '/'.join(d)
    path = os.getcwd()
    if not d in os.getcwd(): os.chdir(d)

    files = [ file for file in os.listdir(os.getcwd()) ]
    write_mode = 'a' if(excel_name in files) else 'w'
    grupos = ret_groups(m, GRUPOS)
    f_totales = Lexico.Frecuencias(DATOS)

    G   = []
    FR  = []
    PAL = []
    IDLs= [ [] for i in range(len(grupos)) ]

    f_i = []
    f = []
    
    for g in grupos:
        G.append(g)
        aux = Lexico.Frecuencias(DATOS, g)
        tmp_n = []
        tmp_p = []
        for fn in aux:
            f.append( fn )
            tmp_n.append( fn[0] )
            tmp_p.append( fn[1] )
        PAL.append(tmp_p)
        FR.append(tmp_n)
        f_i.append(aux)
    f.sort(reverse=True, key=lambda val: val[0] )
    for i in range(15):
        vocablo = f[i][1]
        for g in range(len(grupos)):
            idl = Lexico.IDL_i(vocablo, f_i[g], len(indices))
            IDLs[g].append(idl)

    df = pd.DataFrame()
    for i in range(len(grupos)):
        if len(grupos[i]) < args.max_group: grupos[i].append(-1)
        df["Grupo "+str(i)] = grupos[i]
    df["Valor F"] = [F for i in range(args.max_group)]

    writer = pd.ExcelWriter(excel_name, engine='openpyxl', mode=write_mode)
    df.to_excel(writer, sheet_name='Grupos({})'.format(iteraciones) )
    writer.save()

    df = pd.DataFrame()
    aux = [ [] for i in range(len(grupos)) ]
    for i in range(len(grupos)):
        df["Palabras ({})".format(i)] = PAL[i][:15]
        df["Frecuencias ({})".format(i)] = FR[i][:15]
        for j in range(len(IDLs[i])):
            aux[i].append(IDLs[i][j])

    df["P. más frecuentes"] = [ f[i][1] for i in range(15) ]
    df["Mejores frecuentas"] = [ f[i][0] for i in range(15) ]
    t = []
    for g in range(len(grupos)):
        df["IDL grupo "+str(g)] = aux[g]
    
    for i in range(15):
        v = f[i][1]
        t.append( Lexico.IDL_i(v, f_totales, len(m)) )

    df["Frecuencias totales"] = t

    writer = pd.ExcelWriter(excel_name, engine='openpyxl', mode='a')
    df.to_excel(writer, sheet_name='Palabras({})'.format(iteraciones) )
    writer.save()
    os.chdir(path)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Usar Bacterias :D')
    parser.add_argument('--elit', dest='elitista', action='store_true')
    parser.add_argument('--no-elit', dest='elitista', action='store_false')
    parser.set_defaults(elitista=True)
    parser.add_argument("-F", "--fitness", type=str, choices=["normal", "lexico"], default="normal", help="Tipo de fitness para utilizar. (default = fitness)")
    parser.add_argument("-path", type=str, default="./data/notas.pkl")
    parser.add_argument("-its", type=int, default=2000, help="Número de iteraciones")
    parser.add_argument("-ps","--pob_size", type=int, default=100, help="Tamaño de la población (default = 100) ")
    parser.add_argument("-chp","--change_pob", type=bool, default=False, help="Guardar cambiando población (default = False) ")
    parser.add_argument("-c", "--PC", type=float, default=0.6, help="Probabilidad de Conjugación (default = 0.6) ")
    parser.add_argument("-cl", "--PCL", type=float, default=0.25, help="Probabilidad de Conjugación Libre (default = 0.25) ")
    parser.add_argument('-t', '--PT', type=float, default=0.6, help="Probabilidad de Transformación (default = 0.6) ")
    parser.add_argument('-tl', '--PTL', type=float, default=0.25, help="Probabilidad de Transformación Libre  (default = 0.25) ")
    parser.add_argument('-m', '--PM', type=float, default=0.6, help="Probabilidad de Mutación  (default = 0.6) ")
    parser.add_argument('-ml', '--PML', type=float, default=0.25, help="Probabilidad de Mutación Libre  (default = 0.25) ")
    parser.add_argument('-pl', '--PL', type=float, default=0.5, help="Probabilidad de liberar material genético en el ambiente  (default = 0.5) ")
    parser.add_argument('-mg', '--max_group', type=int, default=5, choices=[3, 4, 5, 6, 7], help="Cantidad de alumnos que tiene el grupo más grande (default = 5) ")
    parser.add_argument('-s', '--stuck', type=bool, default=False, choices=[True, False], help="Indica si se disminuye la población en base al estancamiento  (default = False) ")
    parser.add_argument('-ik', '--iterations_to_kill', type=int, default=50, help="Número de iteraciones que deben pasar para disminuir la población. \n \
                        Si '-s' está activado, representa la cantidad de iteraciones que la población debe estar estancada. (default = 50) ")    
    parser.add_argument('-k', '--kill', type=str, default='f', choices=['b', 'f'], help="Criterio para dejar población viva.\n \
                        \t b (best): Deja viva solo a la mejor bacteria. \n \
                        \t f (fitness): Deja viva a una sola bactereria por fitness. (default)")
    parser.add_argument('-np', '--num_prueba', type=int, default=1, help="Número de la prueba (default = 1)")
    parser.add_argument('-si', '--save_idl', type=str, default='IDL.xlsx', help='Nombre del archivo donde se guardarán los IDL.')
    parser.add_argument('-sf', '--save_fit', type=str, default='BE_lexico.xlsx')
    args = parser.parse_args()

    """
    DEFAULT:
    - elitista = True
    - fitness = "normal" (paper AG)
    - path = "./data/notas.pkl"
    - pob_size = 100
    - change_pob = False
    - PC = 0.6 ; PCL = 0.25
    - PT = 0.6 ; PTL = 0.25
    - PM = 0.6 ; PML = 0.25
    - PL = 0.5
    - max_group = 5
    - stuck = False
    - iterations_to_kill = 50
    - kill = 'f' (Deja una bacteria viva por cada fitness)
    - num_prueba = 1
    - save_idl = "IDL.xlsx"
    - save_fit = "BE_lexico.xlsx"
    """

    # ----- Lectura de datos ----- #
    DATOS = pickle.load( open(args.path, 'rb'), encoding='bytes' )
    indices = [ i for i in DATOS.keys() ]
    _, GRUPOS= crear_grupos(indices, args.max_group) # 4 -> tamaño del grupo más grande
    ALUMN = len(DATOS)
    EVALS = [7, 4, 3, 7] # solo usados para comparar con AG

    ITS = [2000]    
    REPETICIONES = 15


    mejores = UsarBacteria(EVALS, DATOS, GRUPOS, ALUMN, REPETICIONES, args.its, args, True)
    if args.fitness == "lexico": Guardar(args.save_idl, mejores[0], mejores[1], DATOS, args.its)