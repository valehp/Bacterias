import numpy as np
import random
from itertools import compress 
import time

# ======================================================= #
# ==================== Bacteria base ==================== #
# ======================================================= #

class BactariaBase:
	def __init__(self, evaluaciones, notas, grupos, total, psize, PC, PCL, PT, PTL, PM, PML, PL, fitness, max_group=5, kill_percent=-1):
		self.evaluaciones = evaluaciones
		self.NOTAS = notas                  # matriz de notas
		self.indices_alumnos = []
		for i in notas.keys():
		    self.indices_alumnos.append(i)
		self.indices_grupos = grupos        # índices que indican donde empieza cada grupo
		self.total_alumnos = total
		self.custom_fitness = fitness
		self.N = psize                      # Size de la población
		self.pc = PC                        # Probabilidad de conjugación
		self.pc_l = PCL                     # prob conjugacion libre
		self.pt = PT                        # Probabilidad de transformación
		self.pt_l = PTL                     # prob tansformacion libre
		self.pm = PM                        # Probabilidad de mutación
		self.pm_l = PML                     # prob mutación libe
		self.prob_liberar = PL              # Probabilidad de liberar material genético en el ambiente al morir
		self.kill_p = kill_percent
		self.GenerarMinimos(max_group)
		self.reset()
        
	def reset(self):
		self.resistentes = [False for i in range(self.N)]
		self.no_resist = [False for i in range(self.N)]
		self.indices = list(range(self.N))
		self.F = [0] * self.N               # Fitness de la población
		self.B_A = []                       # Aquí se guarda la bacteria que representa el antibiotico
		self.F_A = 0                        # Fitness del antibiótico para no calcularlo de nuevo
		self.P = []                         # Aquí se guarda la población
		self.muertas = []                   # Aquí se guardan los pedazos de material genéticos que se liberan al morir
		self.best_bacteria = []             # Aquí se guarda la mejor solución hasta el momento
		self.best_fitness = float('inf')    # Aquí se guarda el fitness de best_bacteria
		self.max_fitness = -1               # Máximo fitness (o el peor) para comparar si se estanco o no
		#self.MINF = []                      # Aquí se guardan los minimos fitness para el grupo con size i
		self.debug = False
		self.memory = {}


# -------------------- Return de soluciones -------------------- #
	def get_BestSolution(self): # Devuelve la mejor solución entre lo que hay en la población
		f = min(self.F)
		i = self.F.index(f)
		return (self.P[i], f)

	def get_BestAll(self): # Retorna el self.best_fitness y self.best_bacteria
		return (self.best_bacteria, self.best_fitness)

# ----------------------- Inicialización ----------------------- #
	def generar_poblacion(self, N=None):
		if not N: N = self.N
		nums = self.indices_alumnos.copy()
		pobla = []
		for i in range(N):
			aux = random.sample(nums, len(nums))
			pobla.append(aux)
		for bac in pobla:
			self.P.append(bac)


	def GenerarMinimos(self, max_group):
		self.MINF = [ [] for i in range(max_group+1)]
		for GSIZE in range(3, len(self.MINF)):
			aux = []
			for i in range(len(self.evaluaciones)):
				dif = abs(self.evaluaciones[i] - GSIZE)
				if self.evaluaciones[i] >= GSIZE: tmp = GSIZE
				else:
					div = GSIZE // self.evaluaciones[i]
					dif = GSIZE - self.evaluaciones[i]*div
					tmp = (div**3) * (self.evaluaciones[i] - dif) + ( (div+1)**3 )*dif
				aux.append(tmp)
			self.MINF[GSIZE] = aux

# -------------------- -------------------- -------------------- #
	def test_sum(self, T, L):
		# L: Lista con len = n° de opciones (o res posibles de evaluacion i), con todos = 0
		# T: lista de notas del test i
		for i in T: L[i-1] += 1 # Llena L con la cantidad de alumnos con nota L[i]
		for i in range(len(L)): L[i] = L[i] ** 3 # Eleva cada res a 3 (cosa de fitness)
		return sum(L)


	def fitness_group(self, group): # Devuelve la suma de el fitness por grupo
		# group: 1 grupo específico de alumnos
		suma_notas = []
		for i in range(len(self.evaluaciones)):
			aux = []
			for e in group:
				aux.append(self.NOTAS[e][i])
			suma_notas.append(aux)
		aux = 0
		min_ti = self.MINF[len(group)]
		for i in range(len(self.evaluaciones)):
			suma_notas[i] = self.test_sum(suma_notas[i], [0] * self.evaluaciones[i] )
			aux += abs(suma_notas[i] - min_ti[i])
		return aux


    # Métodos para calcular el fitness de un individuo      
	def global_fitness(self): 				# Recorre toda la poblacion y manda a calcular todos los fitness
		self.F = []
		self.indices = list(range(len(self.P)))
		for bacteria in self.P:
			f = self.Fitness(bacteria)
			if f < self.best_fitness:
				self.best_fitness = f
				self.best_bacteria = bacteria
			self.F.append(f)
            

	def Fitness(self, bacteria):   			# Calcula el fitness por bacteria -> separa la bacteria en grupos
		fgrupitos = 0
		all_groups = []
		for index in range(1, len(self.indices_grupos)+1):
			if index >= len(self.indices_grupos):
				i = self.indices_grupos[index-1]
				j = len(bacteria)
			else: 
				i = self.indices_grupos[index-1]
				j = self.indices_grupos[index]
			group = bacteria[i:j]
			if self.custom_fitness: all_groups.append(group)
			else:
				FG = self.fitness_group(group)
				fgrupitos += FG
		if self.custom_fitness:
			FG, self.memory = self.custom_fitness(all_groups, self.NOTAS, self.memory)
			fgrupitos += FG
		return fgrupitos

	def check_fitness(self, it=-1):
		if len(self.F) == 0 or len(self.P) == 0:
			self.generar_poblacion()
			self.global_fitness()
		minimo = min(self.F)
		maximo = max(self.F)
		estancado = False
		if minimo == maximo: estancado = True
		cont = 0
		while minimo == maximo:
			cont += 1
			self.VGL()
			minimo = min(self.F)
			maximo = max(self.F)
			if cont >= 5: break
		return estancado


# -------------------- -------------------- -------------------- #
	# 	VARIACIÓN GENÉTICA 	#
	def _Conjugacion(self, receptora, b_conj, libre):
		"""
		- receptora: bacteria a conjugar
		- b_conj: listado de índices de bacterias disponibles para conjugar
		"""	
		index = random.choice(b_conj)
		donadora = self.P[index]
		nueva = self.combinar_bacterias(donadora, receptora)
		f     = self.Fitness(nueva)

		if libre:
			nueva_l = self.combinar_bacterias(receptora, donadora)
			f_l     = self.Fitness(nueva_l)
			return (nueva, f), (nueva_l, f_l), index

		return nueva, f, index


	def combinar_bacterias(self, donadora, receptora):
		if len(donadora) != len(receptora):
			i = random.randint(0, abs(len(receptora)-len(donadora)-1) )
			j = i+len(donadora)
			seg_don = donadora
		else:
			i = random.randint(0, len(receptora)-3)
			j = random.randint(i+1, len(receptora)-1)
			if j+1 >= len(receptora): seg_don = donadora[i:len(receptora)]
			else: seg_don = donadora[i:j+1]
		if j+1 != len(receptora): h = j+1
		else: h = 0
		nueva = []
		numeracion = []
		listo = False
		for _ in range(len(receptora)):
			if h == len(receptora): h = 0
			if receptora[h] not in seg_don: numeracion.append(receptora[h]) 
			h += 1
		h = 0 # índice para crear la "nueva" bacteria -> máx = len(donadora)
		k = 0 # índice del seg_don (segmento donador)
		for aux in range(len(receptora)):
			if h >= i and k < len(seg_don):
				nueva.append(seg_don[k])
				k += 1
			else:
				nueva.append(numeracion[h])
				h += 1
		return nueva


	def aux_mutacion(self, bac, m1, l, ir):
		"""
		bac: bacteria a mutar
		m1 : primer índice del trozo de mutazión (índice de inicio)
		l  : largo del trozo de mutación
		ir : índice aleatorio en donde se incertará el material
		"""
		if m1 == len(bac) or m1+l >= len(bac): m1 = abs(m1 - l)
		m2 = m1 + l
		BM = [] # Bacteria mutada
		if ir == 0: BM = bac[m1:m2] + bac[:m1] + bac[m2:]
		elif ir == len(bac)-1: BM = bac[:m1] + bac[m2:] + bac[m1:m2]
		else:
			if ir > m1 and ir>= m2: BM = bac[:m1] + bac[m2:ir] + bac[m1:m2] + bac[ir:]
			elif ir > m1 and ir < m2: BM = bac[:m1] + bac[m2:m2+(ir-m1)] + bac[m1:m2] + bac[m2+(ir-m1):]
			else: BM = bac[:ir] + bac[m1:m2] + bac[ir:m1] + bac[m2:]
		return BM


	def VGL(self): # Variación genética libre -> llamar a métodos de variación genética con libre = True
		self.Conjugacion(True)
		self.Transformacion(True)
		self.Mutacion(True)
		self.global_fitness()


	def Regen(self):
		if len(self.P) == 0:
			self.generar_poblacion()
			return
		if len(self.P) >= self.N: return
		dif = self.N - len(self.P)
		for i in range(dif):
			new = random.choice(self.P).copy()
			x1 = random.randint(0, len(new)-1)
			y1 = random.randint(0, len(new)-1)
			x2 = random.randint(0, len(new)-1)
			y2 = random.randint(0, len(new)-1)
			while True:
				if y1 == x1: y1 = random.randint(0, len(new)-1)
				if y2 == x2: y2 = random.randint(0, len(new)-1)
				if y1 != x1 and y2 != x2: break
			aux = new[x1]
			new[x1] = new[y1]
			new[y1] = aux
			aux = new[x2]
			new[x2] = new[y2]
			new[y2] = aux
			self.P.append(new)
			self.F.append(self.Fitness(new))


# -------------------- -------------------- -------------------- #
	# 	ANTIBIÓTICO 	#

	def crear_antib(self): # Torneo
		i1 = random.randint(0, len(self.P)-1)
		i2 = random.randint(0, len(self.P)-1)
		j1 = random.randint(0, len(self.P)-1)
		j2 = random.randint(0, len(self.P)-1)
		mejor1 = i1 if(self.F[i1] > self.F[i2]) else i2
		mejor2 = j1 if(self.F[j1] > self.F[j2]) else j2
		antibiotico = mejor1 if(self.F[mejor1] < self.F[mejor2]) else mejor2
		self.B_A = self.P[antibiotico]
		self.F_A = self.F[antibiotico]
            
	def clasification(self):
		self.resistentes = [False] * len(self.P)
		self.no_resist = [False] * len(self.P)
		for i in range(len(self.P)):
			if self.F[i] <= self.F_A: self.resistentes[i] = True
			else: self.no_resist[i] = True

	def Aplicar_antibiotico(self):
		self.indices = list(range(len(self.P)))
		eliminar = list(compress(self.indices, self.no_resist))
		vivas = list(compress(self.indices, self.resistentes))
		next_generation = []
		next_f = []
		for i in eliminar:
			if random.random() <= self.prob_liberar:
				i1 = random.randint(0, 17)
				i2 = random.randint(i1, 18)
				material_liberado = self.P[i][i1:i2+1]
				if material_liberado: self.muertas.append(material_liberado)
		for i in vivas:
			next_generation.append(self.P[i])
			next_f.append(self.F[i])
		self.P = next_generation
		self.F = next_f
		self.resistentes = [False for i in range(len(self.P))]
		self.no_resist = [False for i in range(len(self.P))]


# ====================== END BACTERIA BASE ====================== #



# =============================================================== #
# ==================== Bacteria con elitismo ==================== #
# =============================================================== #

class BacteriaElitista(BactariaBase):
	def __init__(self, evaluaciones, notas, grupos, total, psize, PC, PCL, PT, PTL, PM, PML, PL, fitness, max_group=5, kill_percent=-1):
		BactariaBase.__init__(self, evaluaciones, notas, grupos, total, psize, PC, PCL, PT, PTL, PM, PML, PL, fitness, max_group, kill_percent)


# --------------- Métodos de verificación y disminución de la población --------------- #
	def kill(self, distint=False):
	# kill se usa para matar p% de la población, SOLO CUANDO ESTÁ TOTALMENTE ESTANCADA => todas las bacterias tienen mismo fitness
	# distint -> True:  se guarda una única bacteria por cada fitnes distinto y las demás son eliminadas.
	# distint -> False: se eliminan todas las bacterias excepto por la que tenga mejor fitness en ese momento.
		if self.kill_p != -1:
			vivas = self.N - int(self.N*self.kill_p)
			if vivas <= 0: vivas = 1
			self.P = self.P[:vivas]
		elif distint:
			if len(self.P) < len(self.F): self.P = self.P[:len(self.F)]
			fitness = np.unique(self.F).tolist()
			vivas = len(fitness)
			p = []
			idx = 0
			for f in fitness:
				idx = self.F.index(f)
				p.append(self.P[idx])
			self.P.clear()
			self.P = p
		elif not distint:
			vivas = 1
			alive = self.get_BestSolution()
			self.P = []
			self.P.append(alive[0])
		dif = self.N - vivas
		self.generar_poblacion(dif)


	def population_check(self, mult):
		if len(self.P) < self.N*mult: return False
		self.kill()        
		return True


# ------------------------ Métodos de variación genética ------------------------- #
	def Conjugacion(self, Libre=False):
		PROB 		  = self.pc_l if(Libre) else self.pc
		bacterias 	  = list(range(len(self.P))) if(Libre) else list(compress(self.indices, self.no_resist))
		para_conjugar = list(range(len(self.P))) if(Libre) else list(compress(self.indices, self.resistentes))
		if not para_conjugar: return

		for index in bacterias:
			b = self.P[index]
			prob = random.random()
			if prob <= PROB:
				if Libre:
					(nb1, f1), (nb2, f2), i = self._Conjugacion(b, para_conjugar, Libre)
					self.P[index] = nb1
					self.F[index] = f1
					self.P[i] = nb2
					self.F[i] = f2

				else:
					nb, f, i = self._Conjugacion(b, para_conjugar, Libre)
					self.P[index] = nb
					self.F[index] = f
				"""
				if Libre:
					(new1, _), (new2, _), _ = self._Conjugacion(b, para_conjugar, Libre)
					self.P.append(new1)
					self.P.append(new2)
				else:
					new, _, _ = self._Conjugacion(b, para_conjugar, Libre)
					self.P.append(new)
				"""


	def Transformacion(self, Libre=False):
		if not self.muertas: return
		if Libre:
			PROB = self.pt_l
			bacterias = list(range(len(self.P)))
		else:
			PROB = self.pt
			bacterias = list(compress(self.indices, self.no_resist))
		for index in bacterias:
			if not self.muertas: break
			b = self.P[index]
			prob = random.random()
			if prob <= PROB:
				material = random.choice(self.muertas)
				nueva = self.combinar_bacterias(material, b)
				nueva_fit = self.Fitness(nueva)
				self.P[index] = nueva
				self.F[index] = nueva_fit
				self.muertas.remove(material)
				if not Libre:
					if nueva_fit < self.F_A:
						self.no_resist[index] = False
						self.resistentes[index] = True
			"""
			if prob <= PROB:
				material = random.choice(self.muertas)
				nueva = self.combinar_bacterias(material, b)
				self.P.append(nueva)
				self.muertas.remove(material)
			"""


	def Mutacion(self, Libre=False):
		if Libre:
			PROB = self.pm_l
			indx = list(range(len(self.P)))
		else:
			PROB = self.pm
			indx = list(range(self.N))
		#print("Mutacion -> inx: ", len(indx), " - restist: ", len(self.resistentes), " - no resist: ", len(self.no_resist), " - bacs: ", len(self.P))
		for i in indx:
			prob = random.random()
			if prob <= PROB:

				mutada = self.P[i].copy()
				m1 = random.randint(0, len(mutada)-1)
				l = random.randint(1, len(mutada)-2)
				ir = random.randint(0, len(mutada)-1) # indice random donde incertar la mutación
				while ir == m1: ir = random.randint(0, len(mutada)-1)
				nueva_bac = self.aux_mutacion(mutada, m1, l, ir)
				nueva_fit = self.Fitness(nueva_bac)
				if Libre:
					self.P[i] = nueva_bac
					self.F[i] = nueva_fit
				else:
					if self.no_resist[i]:
						self.P[i] = nueva_bac
						self.F[i] = nueva_fit
						if nueva_fit < self.F_A:
							self.no_resist[i] = False
							self.resistentes[i] = True
					elif self.resistentes[i] and self.F[i] > nueva_fit:
						self.P[i] = nueva_bac
						self.F[i] = nueva_fit
				#if Libre or self.no_resist[i]: self.P.append(nueva_bac)
				#elif self.resistentes[i] and nueva_fit < self.F[i]: self.P.append(nueva_bac)



	def use(self, iterations=1000, kill_iters=50, distint=True):
		fits = []
		times = []
		start = time.time()
		i, last_k = 0, 0

		self.reset()
		self.generar_poblacion()

		for it in range(iterations):
			if last_k - i == kill_iters:
				last_k = i
				self.kill(distint)

			self.global_fitness()
			self.crear_antib()
			self.clasification()

			self.Conjugacion()
			self.Transformacion()
			self.Mutacion()

			self.global_fitness()
			self.clasification()
			self.Aplicar_antibiotico()
			self.Regen()
			self.check_fitness()

			fits.append( self.get_BestSolution()[1] )
			times.append( time.time() - start )

		return fits, times, self.get_BestSolution()



# -------------------- END BACTERIA ELITISTA -------------------- #





# =============================================================== #
# ==================== Bacteria sin elitismo ==================== #
# =============================================================== #

class Bacteria(BactariaBase):
	def __init__(self, evaluaciones, notas, grupos, total, psize, PC, PCL, PT, PTL, PM, PML, PL, fitness, max_group=5, kill_percent=-1):
		BactariaBase.__init__(self, evaluaciones, notas, grupos, total, psize, PC, PCL, PT, PTL, PM, PML, PL, fitness, max_group, kill_percent)



# ------------------------ Métodos de variación genética ------------------------- #
	def Conjugacion(self, Libre=False):
		PROB 		  = self.pc_l if(Libre) else self.pc
		bacterias 	  = list(range(len(self.P))) if(Libre) else list(compress(self.indices, self.no_resist))
		para_conjugar = list(range(len(self.P))) if(Libre) else list(compress(self.indices, self.resistentes))
		if not para_conjugar: return

		for index in bacterias:
			b = self.P[index]
			prob = random.random()
			if prob <= PROB:
				if Libre:
					(nb1, f1), (nb2, f2), i = self._Conjugacion(b, para_conjugar, Libre)
					self.P[index] = nb1
					self.F[index] = f1
					self.P[i] = nb2
					self.F[i] = f2

				else:
					nb, f, i = self._Conjugacion(b, para_conjugar, Libre)
					self.P[index] = nb
					self.F[index] = f

    
	def Transformacion(self, Libre = False):
		if not self.muertas: return
		if Libre:
			PROB = self.pt_l
			bacterias = self.indices
		else:
			PROB = self.pt
			bacterias = list(compress(self.indices, self.no_resist))
		for index in bacterias:
			if not self.muertas: break
			b = self.P[index]
			prob = random.random()
			if prob <= PROB:
				material = random.choice(self.muertas)
				nueva = self.combinar_bacterias(material, b)
				nueva_fit = self.Fitness(nueva)
				self.P[index] = nueva
				self.F[index] = nueva_fit
				self.muertas.remove(material)
				if not Libre:
					if nueva_fit < self.F_A:
						self.no_resist[index] = False
						self.resistentes[index] = True
    
	def Mutacion(self, Libre = False):
		if Libre: PROB = self.pm_l
		else: PROB = self.pm
		for i in range(self.N):
			prob = random.random()
			if prob <= PROB:
				mutada = self.P[i].copy()
				m1 = random.randint(0, len(mutada)-1)
				l  = random.randint(1, len(mutada)-2)
				ir = random.randint(0, len(mutada)-1) # indice random donde incertar la mutación
				while ir == m1: ir = random.randint(0, len(mutada)-1)
				nueva_bac = self.aux_mutacion(mutada, m1, l, ir)
				nueva_fit = self.Fitness(nueva_bac)
				if Libre:
					self.P[i] = nueva_bac
					self.F[i] = nueva_fit
				else:
					if self.no_resist[i]:
						self.P[i] = nueva_bac
						self.F[i] = nueva_fit
						if nueva_fit < self.F_A:
							self.no_resist[i] = False
							self.resistentes[i] = True
					elif self.resistentes[i] and self.F[i] > nueva_fit:
						self.P[i] = nueva_bac
						self.F[i] = nueva_fit



	def use(self, iterations=1000):
		fits = []
		times = []
		start = time.time()

		self.reset()
		self.generar_poblacion()

		for it in range(iterations):
			self.global_fitness()
			self.crear_antib()
			self.clasification()

			self.Conjugacion()
			self.Transformacion()
			self.Mutacion()

			self.global_fitness()
			self.clasification()
			self.Aplicar_antibiotico()
			self.Regen()
			self.check_fitness()

			fits.append( self.get_BestSolution()[1] )
			times.append( time.time() - start )

		return fits, times, self.get_BestSolution()