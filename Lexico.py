# -*- coding: utf-8 -*-
import numpy as np 
import pickle
import time
from itertools import permutations

def distance(str1, str2):
	d = dict()
	for i in range(len(str1)+1):
		d[i] = dict()
		d[i][0] = i

	for i in range(len(str2)+1):
		d[0][i] = i
	
	for i in range(1, len(str1)+1):
		for j in range(1, len(str2)+1):
			d[i][j] = min(d[i][j-1]+1, d[i-1][j]+1, d[i-1][j-1]+(not str1[i-1] == str2[j-1]))
	return d[len(str1)][len(str2)]


def HET(grupo, pals, d): # (1)
	suma = 0
	for i in range( len(grupo) ):
		for j in range( i+1, len(grupo) ):
			p1, p2 = grupo[i], grupo[j]
			tupla  = (p1, p2)
			if tupla not in d: d[tupla] = []
			for l in range(len(pals[p1])):
				if tupla in d and len(d[tupla]) > l: suma += d[tupla][l]
				else:
					#aux = levenshtein( pals[p1][l], pals[p2][l] )
					aux = distance( pals[p1][l], pals[p2][l]  )
					suma += aux
					d[tupla].append(aux)
	return (suma, d)


def HOM(grupos, pals, d): # (2)
	suma = 0
	for i in range( len(grupos) ):
		for j in range( i+1, len(grupos) ):
			n, m = tuple(grupos[i]), tuple(grupos[j])
			
			if n not in d:
				H_n, d = HET(n, pals, d)
				d[n]   = H_n
			else: H_n = d[n]

			if m not in d: 
				H_m, d = HET(m, pals, d)
				d[m]   = H_m
			else: H_m = d[m] 
			
			suma += abs( H_n - H_m )
	return (suma, d)


def F(grupos, lex, d={}): # (3)
	suma = 0
	hom, d = HOM(grupos, lex, d)
	for grupo in grupos:
		aux = d[ tuple(grupo) ] 
		suma  += aux
	total= (1/hom) * suma
	return (total, d)


def CreateGroups(g, max_len):
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


def Frecuencias(lexico, grupos=[]):
	fn = []
	f  = { 'NUM':{}, 'CALC':{}, 'ESTR':{}, 'GEOM':{}, 'AZAR':{} }
	if not grupos:
		for key, values in lexico.items():
			for l in range(len(values)):
				if l == 0: p = 'NUM'
				if l == 1: p = 'CALC'
				if l == 2: p = 'ESTR'
				if l == 3: p = 'GEOM'
				if l == 4: p = 'AZAR'
				for palabra in values[l]:
					if palabra not in f[p]: f[p][palabra] = 1
					else: f[p][palabra] += 1
	else:
		for indice in grupos:
			values = lexico[indice]
			for l in range(len(values)):
				if l == 0: p = 'NUM'
				if l == 1: p = 'CALC'
				if l == 2: p = 'ESTR'
				if l == 3: p = 'GEOM'
				if l == 4: p = 'AZAR'
				for palabra in values[l]:
					if palabra not in f[p]: f[p][palabra] = 1
					else: f[p][palabra] += 1

	for k, dicts in f.items():
		for key, frec in dicts.items():
			fn.append( (frec, key) )
	
	fn.sort(reverse=True, key=lambda val: val[0] )
	return fn 


def IDL_vi (vocablo, listas, Lambda=0.9):
    f = []
    ml = 0
    for l in listas: ml = len(l) if (len(l) > ml) else ml
    
    f = [ 0 for i in range(ml) ]
    
    for l in listas:
        for i in range(len(l)):
            p = l[i]
            if p == vocablo: f[i] += 1
    
    #print("Frecuencias para ", vocablo, "\n", f)
    
    idl = 0
    for i in range(len(f)):
        idl += f[i] * (Lambda ** i)
    idl /= len(listas)
    return idl
    

def IDL_i(v_i, f, N, Lambda=0.9):
	"""
	- fn: Listado de frecuencias de las palabras
	- N : total de individuos
	- lambda: factor lambda de ponderación
	"""
	i = 0
	fn = 0
	f_ant = -1
	for j in range(len(f)):
		if f_ant == -1: f_ant = f[j][0]
		elif f_ant != f[j][0]: i+= 1

		if f[j][1] == v_i:
			fn = f[j][0]
			break

		f_ant = f[j][0]

	idl = (fn * Lambda**i)/N
	return idl


def IDL(fn, N, Lambda=0.9):
	"""
	- fn: Listado de frecuencias de las palabras
	- N : total de individuos
	- lambda: factor lambda de ponderación
	"""
	idl = 0.0
	print( len(fn) )
	for i in range(len(fn)):
		idl += ( fn[i][0] * Lambda**i )
	idl /= N
	return idl


class AUTO:
	def __init__(self, alumnos, lexico, grupos):
		"""
		- alumnos: Lista con índices correspondientes a los alumnos
		- lexico : Diccionario con grupos[indice] = listas con el vocabulario de los alumnos
		- grupos : índices de los grupos 
		"""
		self.alumnos = alumnos 
		self.lexico  = lexico
		self.indices_grupos = grupos
		self.best_F  = float('inf')
		self.best 	 = []
		self.idl 	 = []

	def GetGrupos(self, lista):
		all_groups = []
		for index in range(1, len(self.indices_grupos)+1):
			if index >= len(self.indices_grupos):
				i = self.indices_grupos[index-1]
				j = len(lista)
			else: 
				i = self.indices_grupos[index-1]
				j = self.indices_grupos[index]
			group = lista[i:j]
			all_groups.append(group)
		return all_groups

	def Buscar(self, lista):
		aux  = {}
		cont = 0

		for p in permutations(self.alumnos):
			if cont % 100 == 0: print(cont, end=" -> ")
			g = self.GetGrupos(p)
			F_p = F(g, self.lexico, aux)
			if F_p < self.best_F:
				self.best = [p]
				self.best_F = F_p
			elif F_p == self.best_F: self.best.append(p)
			cont += 1

		if len(p) > 1 and type(p[0]) == list and type(p[-1]) == list:
			for i in p: idl.append( IDL(i) )
		else: idl.append( IDL(p) )
		self.idl = idl
		return (self.best, self.best_F, idl)