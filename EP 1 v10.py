import numpy as np
import matplotlib.pyplot as plt
import time

#Vitor Rana Camarotto      nusp: 103 380 51

# A função main está no final do programa

#===========|| Funções de Soluções exatas ||==================

def u_exataB(x,t):
    """ Calcula a função exata da letra b)"""
    u = np.exp(t-x)*np.cos(5*t*x)
    return u

def u_exataA2(x,t):
    """Calcula a função exata da letra a)"""
    u = (1 + np.sin(10*t))*x*x*(1-x)*(1-x)
    return u

def Us_exata(N,Resposta,M):
    """Função que devolve uma lista 'USS' para plotagem do gráfico da solução exata desejada"""
    t_exato = np.linspace(0,1,M)
    x_exato = np.linspace(0,1,N)
    USS = []
    for i in t_exato:
        Us_exata = []
        for j in x_exato:
            if Resposta == "a":
                uEx = u_exataA2(j,i)
            elif Resposta == "b":
                uEx = u_exataB(j,i)
            Us_exata.append(uEx)
        USS.append(Us_exata)
    return USS

#====|| Funções das condições de fronteira da letra b) ========

def g1b(t):
    "Calcula a condição de fronteira para x = 0, t variando"
    g1 = np.exp(t)
    return g1
def g2b(t):
    "Calcula a condição de fronteira para x = 1 t variando"
    g2 = np.exp(t-1)*np.cos(5*t)
    return g2


#===============|| Funções de fontes f(x,t) ||==================

def f_fonteA2(x,t):
    "Função fonte f(t,x) da letra a)"
    f = 10*np.cos(10*t)*x*x*(1-x)*(1-x) - (1 + np.sin(10*t))*(12*x**2 - 12*x +2)
    return f

def f_fonteB(x,t):
    "Função fonte f(t,x) da letra b)"
    f = 5*np.exp(t-x)*(np.sin(5*t*x)*(-x-2*t)+np.cos(5*t*x)*5*t**2)
    return f
#==========================================================
#==========================================================



# ======|| Funções da primeira tarefa, letra c) ||=============

def r(t):
    "Função r(t)"
    r = 10000*(1-2*t*t)
    return r

def gh(x,p):
    "Função de integração"
    h = dx
    Gh = []
    for i in x:
        if p - (h/2) <= i <= p + (h/2):
            g = 1/h
            Gh.append(g)
        else:
            g = 0
            Gh.append(g)
    return Gh

def f_fonteC(x,k,p):
    "Função fonte f(t,x) da letra c)"
    Gh = gh(x,p)
    F =[]
    for i in Gh:
        f = r(k) * i
        F.append(f)
    return F

#==========================================================
#==========================================================

#=====|| Função que define condições iniciais ||===========

def cond_iniciais(n,u_a,x):  # n identifica a tarefa; u_a é uma lista que representa o eixo x
    "Função que devolve lista com as condições iniciais"
    if n == 2: #Letra B
        if prime == 1:
            for i in x:
                k = int(i*(N-1))
                u_a[k] = u_exataB(i,0)
        else:
            for k in range(len(x)):
                u_a[k] = u_exataB(x[k],0)
    elif n == 1: #Letra A2
        if prime == 1:
            for i in x:
                k = int(i*(N-1))
                u_a[k] = (i**2)*(1-i)**2
        else:
            for k in range(len(x)):
                u_a[k] = u_exataA2(x[k],0)
    return u_a
#==========================================================
#==========================================================


#========|| Função principal da primeira tarefa ||=========

def Dif_finita(x,t,u_a,u,dx,dt,letra):
    """Função principal da primeira tarefa"""
    """Recebe duas listas com dispersões em x e em t; e listas que guardam os valores das temperaturas u_a (atual) e u (próximo t)"""
    global Passos
    Passos = 1
    g = letra
    X = x
    Us = []
    u = np.array(u)
    u_a = np.array(u_a)
    if g == "a":
        """resolve letra a) pelo metodo de diferencas finitas"""
        u_a = cond_iniciais(1,u_a,x)
        for k in range(0,len(t)):
            Passos += 1
            u[1:N-1] = u_a[1:N-1] + dt*( (u_a[0:N-2] - 2*u_a[1:N-1] + u_a[2:N]) /(dx**2) + f_fonteA2(x[1:N-1],t[k]))
            u[0] = 0
            u[-1] = 0 
            Us.append(u)
            u = np.round(u,12)
            u_a, u = u, u_a

    elif g == "b":
        """resolve letra b) pelo metodo de diferencas finitas"""
        u_a = cond_iniciais(2,u_a,X)

        for k in range(0,len(t)):
            Passos += 1
            u[1:N-1] = u_a[1:N-1] + dt*( (u_a[0:N-2] - 2*u_a[1:N-1] + u_a[2:N]) /(dx**2) + f_fonteB(x[1:N-1],t[k]))
            u[0] = np.exp(t[k])
            u[-1] = np.exp(t[k]-1)*np.cos(5*t[k])
            Us.append(u)
            u = np.round(u,12)
            u_a, u = u, u_a

    elif g == "c":
        """resolve letra c) pelo metodo de diferencas finitas"""
        pe = 0.25
        for k in range(0,len(t)):
            Passos += 1
            u[1:N-1] = u_a[1:N-1] + dt*( (u_a[0:N-2] - 2*u_a[1:N-1] + u_a[2:N]) /(dx**2) + f_fonteC(x[1:N-1],t[k],pe))
            u[0] = 0
            u[-1] = 0
            Us.append(u)
            u = np.round(u,12)
            u_a, u = u, u_a
    return Us,u_a

#==========================================================
#==========================================================


#=======|| Funções da Segunda Tarefa ||=====================

def LDLt(A):
    """Decompõe a matriz A em 3 matrizes, tal que A = (L)*(D)*(Lt)
            Devolve as três matrizes L,D,Lt, e os vetores l (constituinte das matrizes L e Lt) e d (constituinte da matriz D)"""

    n = len(A)
    l = np.zeros(n)
    u = np.zeros(n)
    C = np.zeros(n)
    d = np.zeros(n)
    D = np.zeros((n,n))
    L = np.zeros((n,n))
    Lt = np.transpose(L)

    for i in range(n-1):
        C[i] = A[i][i+1]

    u[0] = A[0][0]
    for i in range(1,n):
        l[i] = A[i][i-1] / u[i-1]
        u[i] = A[i][i] - l[i] * C[i-1]

    """Definir matriz L, por motivos de teste"""
    for i in range(n):
        for j in range(n):
            if i == j:
                L[i][i] = 1
            else:
                if i == j + 1:
                    L[i][j] = l[i]

    """Definir matriz D, por motivos de teste"""
    for i in range(n):
        for j in range(n):
            if i == j:
                D[i][i] = u[i]
                d[i] = u[i]

    return L,D,Lt,l,d

def ResolveSist(l,d,u_a,b,t,k):
    """Função que resolve o sistema linear (A)x=b = (LDLt)x = b"""
    """matrizes auxiliares para resolução do sistema"""
    Q = len(b)
    y = np.zeros(Q)
    z = np.zeros(Q)
    X = np.zeros(Q)
    for i in range(0,Q):
        if i == 0:
            y[i] = b[i]
        else:
            y[i] = b[i] - l[i] * y[i-1]

    for i in range(0,Q):
        z[i] = y[i] / d[i]

    for i in range(Q-1,-1,-1):
        if i == Q-1:
            X[i] = z[i]
        else:
            X[i] = z[i] - l[i+1] * X[i+1]

    #Definir condiçoes de fronteira
    if letra == "a":
        X = np.append(X,u_exataA2(0,t[k+1]))
        X = np.insert(X,[0],u_exataA2(1,t[k+1]))
    elif letra == "b":
        X = np.append(X,u_exataB(1,t[k+1]))
        X = np.insert(X,[0],u_exataB(0,t[k+1]))
    elif letra == "c":
        X = np.append(X,0)
        X = np.insert(X,[0],0)

    return X

def imprime(A):
    """imprime uma matriz A"""
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] = np.round(A[i][j],3)
            print(A[i][j], end ='\t')
        print()

#========|| Função principal da segunda tarefa ||=========

def sis_linear(x,t,dx,dt,letra):
    """matriz A dos coeficientes"""
    global Passos
    Passos = 1
    A = np.zeros((N-2,N-2))
    for i in range(0,N-2):
        if i != 0: 
            if metodo == "e":
                A[i,i-1] = -Lambda
            elif metodo == "c":
                A[i,i-1] = -Lambda/2
        if i != N-3:
            if metodo == "e":
                A[i,i+1] = -Lambda
            elif metodo == "c":
                A[i,i+1] = -Lambda/2
        if metodo == "e":
            A[i,i] = 1 + 2*Lambda
        elif metodo == "c":
            A[i,i] = 1 + Lambda

    """Decompor a matriz A"""
    L,D,Lt,l,d = LDLt(A)
    Us = []
    if letra == "a":
        if metodo == "e":
            """Resolve o sistema da letra a) pelo método de Euler"""
            """fatores que compoem b, de b = Ax"""
            x = np.delete(x,0)
            x = np.delete(x,-1)
            u_a = np.zeros(N-2)
            detefes = np.zeros(N-2)
            u_a = cond_iniciais(1,u_a,x)
            for k in range(0,len(t)-1):
                Passos += 1
                detefes = dt*f_fonteA2(x[:],t[k+1])
                b = u_a + detefes
                xx = ResolveSist(l,d,u_a,b,t,k)
                Us.append(xx)
                if k != len(t)-2:
                    xx = np.delete(xx,0)
                    xx = np.delete(xx,-1)
                    u_a, xx = xx, u_a

        elif metodo == "c":
            """Resolve o sistema da letra a) pelo método de Crank-Nicolson"""

            u_a = np.zeros(N)
            u_a = cond_iniciais(1,u_a,x)
            x = np.delete(x,0)
            x = np.delete(x,-1)
            for k in range(0,len(t)-1):
                Passos += 1
                ene = len(u_a)
                detefes = (dt/2)*( f_fonteA2(x[:],t[k]) + f_fonteA2(x[:],t[k+1]) )
                yus = (Lambda/2)*( u_a[0:ene-2] + u_a[2:ene] ) + ((1-Lambda) * u_a[1:ene-1])
                b = np.zeros(N)
                b = detefes + yus
                xx = ResolveSist(l,d,u_a,b,t,k)

                Us.append(xx)
                if k != len(t)-2:
                    u_a, xx = xx, u_a

    elif letra == "b":
        if metodo == "e":
            """Resolve o sistema da letra b) pelo método de Euler"""
            """fatores que compoem b, de b = Ax"""

            x = np.delete(x,0)
            x = np.delete(x,-1)
            u_a = np.zeros(N-2)
            detefes = np.zeros(N-2)
            u_a = cond_iniciais(2,u_a,x)
            lamge1 = 0
            lamge2 = 0
            for k in range(0,len(t)-1):
                Passos += 1
                detefes = dt*f_fonteB(x[:],t[k+1])
                b = u_a + detefes
                lamge1 = Lambda * g1b(t[k+1])
                lamge2 = Lambda * g2b(t[k+1])
                b[0] += lamge1
                b[-1] += lamge2
                xx = ResolveSist(l,d,u_a,b,t,k)
                Us.append(xx)

                if k != len(t)-2:
                    xx = np.delete(xx,0)
                    xx = np.delete(xx,-1)
                    u_a, xx = xx, u_a

        elif metodo == "c":
            """Resolve o sistema da letra b) pelo método de Crank-Nicolson"""

            u_a = np.zeros(N)
            u_a = cond_iniciais(2,u_a,x)
            u_a = np.array(u_a)
            x = np.delete(x,0)
            x = np.delete(x,-1)
            lamge1 = 0
            lamge2 = 0
            for k in range(0,len(t)-1):
                Passos += 1
                ene = len(u_a)
                detefes = (dt/2)*( f_fonteB(x[:],t[k]) + f_fonteB(x[:],t[k+1]) )
                yus = (Lambda/2)*( u_a[0:ene-2] + u_a[2:ene] ) + (1-Lambda) * u_a[1:ene-1]
                b = np.zeros(N)
                b = detefes + yus
                lamge1 = (Lambda/2)*g1b(t[k+1])
                lamge2 = (Lambda/2)*g2b(t[k+1])
                b[0] += lamge1
                b[-1] += lamge2
                xx = ResolveSist(l,d,u_a,b,t,k)
                Us.append(xx)
                if k != len(t)-2:
                    u_a, xx = xx, u_a

    elif letra == "c":
        if metodo == "e":
            """Resolve o sistema da letra b) pelo método de Euler"""
            """fatores que compoem b, de b = Ax"""

            x = np.delete(x,0)
            x = np.delete(x,-1)
            u_a = np.zeros(N-2)
            pe = 0.80
            for k in range(0,len(t)-1):
                Passos += 1
                F = f_fonteC(x[:],t[k+1],pe)
                for i in range(len(F)):
                    F[i] = dt*F[i]
                b = u_a + F
                xx = ResolveSist(l,d,u_a,b,t,k)
                Us.append(xx)
                if k != len(t)-2:
                    xx = np.delete(xx,0)
                    xx = np.delete(xx,-1)
                    u_a, xx = xx, u_a

        elif metodo == "c":
            """Resolve o sistema da letra c) pelo método de Crank-Nicolson"""
            u_a = np.zeros(N)
            u_a = np.array(u_a)
            x = np.delete(x,0)
            x = np.delete(x,-1)
            pe = 0.25
            for k in range(0,len(t)-1):
                Passos += 1
                ene = len(u_a)
                F = f_fonteC(x[:],t[k],pe)
                F1 = f_fonteC(x[:],t[k+1],pe)
                EFE = np.array(F) + np.array(F1)
                for i in range(len(EFE)):
                    EFE[i] = (dt/2)*EFE[i]
                yus = (Lambda/2)*( u_a[0:ene-2] + u_a[2:ene] ) + (1-Lambda) * u_a[1:ene-1]
                yus = np.array(yus)
                b = np.zeros(N)
                b = EFE + yus
                xx = ResolveSist(l,d,u_a,b,t,k)
                Us.append(xx)
                if k != len(t)-2:
                    u_a, xx = xx, u_a

    return Us, xx

def erro(Us, USS):
    """Calcula o erro da aproximação feita"""
    USS = np.array(USS)
    Us = np.array(Us)
    if prime == 2:
        Us = np.delete(Us, (0), axis=0)

    Er = abs(USS[:] - Us[:])

    return Er

#==========================================================
#==================== MAIN FUNCTION =======================
#==========================================================

def main():
    global run_program
    run_program = 1
    while run_program == 1:

        start = time.time()
        global N
        global Lambda
        global dx
        global prime # Identifica a tarefa 1 ou 2
        global letra # Identifica letra do exercicio a ser realizado (a,b,c)
        global metodo # Identifica metodo de Euler ou Crank-Nicolson
        global Metodo # Utilizada para titular gŕaficos com o método utilizado
        global Passos
        prime = int(input("olá, qual tarefa devo realizar? ( 1 / 2 )"))

        if prime == 1:
            qw = 0

            while qw != 1:
                letra = str(input("Qual letra devo realizar? ( a / b / c )"))
                if letra in "abc" :
                    print()
                    print("Certo!")
                    print()
                    print()
                    qw += 1
                else:
                    print("Tente digitar uma das letras indicadas dessa vez")
                    print()
            print("===|| Primeira tarefa - Letra {}) ||===)" .format(letra))
            T = 1
            print("E o invervalo de tempo é T = 1")
            print()
            N = int(input("digite N = "))
            Lambda = float(input("digite o valor de λ = "))
            N += 1
            print()
            print()
            print("Aguarde enquanto os gráficos são calculados")
            x = np.linspace(0,1,N)
            dx = x[1] - x[0]
            dt = Lambda * (dx**2)
            M = int(T/dt)
            t = np.linspace(0,1,M)
            u = np.zeros(N)           #incógnita u a ser calculada no próximo t
            u_a = np.zeros(N)         #valores de u ao longo de x, no t atual
            processa = time.time()
            """ Resolvendo a equação diferencial parcial pelo método das diferenças finitas """
            Us, Sol_final = Dif_finita(x,t,u_a,u,dx,dt,letra)
            # Montando a solução exata
            if letra != "c":
                USS = Us_exata(N,letra,M)

        elif prime == 2:
            qw = 0

            while qw != 1:
                letra = str(input("Qual letra devo realizar? ( a / b / c )"))
                if letra in "abc" :
                    print()
                    print("Certo!")
                    print()
                    print()
                    qw += 1
                else:
                    print("Tente digitar uma das letras indicadas dessa vez")
                    print()
            metodo = input("Qual método de resolução devo utilizar, Euler ou Crank-Nicolson? ( e / c) ")
            print()
            print()
            print()
            print()
            print()
            if metodo == "e":
                print("===|| Segunda Tarefa - Letra {}) - Método de Euler ||===)" .format(letra))
                print()
            elif metodo == "c":
                print("===|| Segunda Tarefa - Letra {}) - Método de Crank-Nicolson ||===)" .format(letra))
                print()
            T = 1
            print("E o invervalo de tempo é T = 1")
            print()
            N = int(input("digite N = "))
            Lambda = N
            N += 1
            print()
            print()
            print("Aguarde enquanto os gráficos são calculados")
            x = np.linspace(0,1,N)
            dx = x[1] - x[0]
            dt = dx
            M = int(T/dt)
            t = np.linspace(0,1,M)
            processa = time.time()
            """ Montando e resolvendo o sistema linear """
            Us, Sol_final = sis_linear(x,t,dx,dt,letra)

        elif prime == 3:
            print("Suas opções são")
            print()
            print("Verificar todos os erros da primeira tarefa (1)")
            print("Verificar a decomposição da matriz A em L D Lt (2)")
            opcao = int(input())
            if opcao == 1: 
                T = 1
                LAMBDAS = [0.25,0.50]
                ENES = [10,20,40,80,160,320]
                LETRA = "ab"
                for l in LAMBDAS:
                    Errito = []
                    for Letra in LETRA:
                        for en in ENES:
                            N = en
                            Lambda = l
                            letra = Letra
                            x = np.linspace(0,1,en)
                            dx = x[1] - x[0]
                            dt = l * (dx**2)
                            M = int(T/dt)
                            t = np.linspace(0,1,M)
                            u = np.zeros(en)
                            u_a = np.zeros(en)

                            Us, Sol_final = Dif_finita(x,t,u_a,u,dx,dt,Letra)
                            USS = Us_exata(N,Letra,M)
                            ERERO = np.amax(erro(USS,Us))
                            Errito.append(ERERO)
                            print("λ = {}, N = {}, Letra {}), o maior erro é {} " .format(l,en,Letra,ERERO))

            elif opcao == 2:
                N = int(input("Digite o valor de N = (Recomendado: 10)"))
                metodo = input("Qual método de resolução devo utilizar, Euler ou Crank-Nicolson? ( e / c) ")
                Lambda = N
                A = np.zeros((N-2,N-2))
                for i in range(0,N-2):
                    if i != 0:
                        if metodo == "e":
                            A[i,i-1] = -Lambda
                        elif metodo == "c":
                            A[i,i-1] = -Lambda/2
                    if i != N-3:
                        if metodo == "e":
                            A[i,i+1] = -Lambda
                        elif metodo == "c":
                            A[i,i+1] = -Lambda/2
                    if metodo == "e":
                        A[i,i] = 1 + 2*Lambda
                    elif metodo == "c":
                        A[i,i] = 1 + Lambda

                """Decompor a matriz A"""
                L,D,Lt,l,d = LDLt(A)
                print()
                print()
                print()
                print("L = ")
                imprime(L)
                print()
                print()
                print("D = ")
                imprime(D)
                print()
                print()
                print("Lt = ")
                imprime(Lt)
                print()

            run_program = 0
            break


        if letra != "c":
                USS = Us_exata(N,letra,M)

#======================================================
#======================================================

#===========|| Plotagem dos gŕáficos ||================

        """ Plotagem dos gráficos a cada 0.1s """
        eixo_t = np.linspace(0,1,N)
        Nepkins = int(M/10)
        if prime == 1:
            for e in range(0,10):
                Fepkins = e*Nepkins
                if e == 0:
                    if letra == "a":
                        u_a = cond_iniciais(1,np.zeros(N),x)
                        plt.plot(eixo_t,u_a,color=(0,0,1,1),label ='T = 0')
                    elif letra == "b":
                        u_a = cond_iniciais(2,np.zeros(N),x)
                        plt.plot(eixo_t,u_a,color=(0.8,0,0.6,0.2), label ='T = 0')
                    elif letra == "c":
                        u_a = np.zeros(N)
                        plt.plot(eixo_t,u_a,color=(0.8,0,0.6,0.2), label ='T = 0')
                else:
                    if letra == "a":
                        if e<5:
                            alpha = e*0.1+0.1
                            plt.plot(eixo_t,Us[Fepkins-1],color = (alpha,0,1,alpha),label='T = {}' .format(e/10))
                        elif e<10:
                            alpha = e*0.1
                            plt.plot(eixo_t,Us[Fepkins-1],color = (0.4,alpha,1-alpha,alpha),label='T = {}' .format(e/10))
                    elif letra == "b":
                        alpha = 0.1*e
                        plt.plot(eixo_t,Us[Fepkins-1],color = (1,0,1,alpha),label='T = {}' .format(e/10))
                    elif letra == "c":
                        alpha = 0.1*e
                        plt.plot(eixo_t,Us[Fepkins-1],color = (0.8,0,0.3,alpha),label='T = {}' .format(e/10))

            if letra == "a":
                        plt.plot(eixo_t,Sol_final,color=(0,1,0,1),label ='T = 1')
            elif letra == "b":
                        plt.plot(eixo_t,Sol_final,color=(0.8,0,0.6,1), label ='T = 1')
            elif letra == "c":
                        plt.plot(eixo_t,Sol_final,color=(0.8,0,0.3,1), label ='T = 1')
            plt.ylabel('Temperatura')
            plt.xlabel('Ponto X da barra')
            plt.title('Gráfico da Tarefa {}, letra {})\nAproximação pelo método das diferenças finitas\n N = {}, λ = {}' .format(prime,letra,N-1,Lambda))
            plt.legend()
            plt.show()
        elif prime == 2:
            for e in range(0,10):
                Fepkins = e*Nepkins
                if e == 0:
                    if letra == "a":
                        u_a = cond_iniciais(1,np.zeros(N),x)
                        plt.plot(eixo_t,u_a,color=(0,0,1,1),label ='T = 0')
                    elif letra == "b":
                        u_a = cond_iniciais(2,np.zeros(N),x)
                        plt.plot(eixo_t,u_a,color=(0.8,0,0.6,0.2), label ='T = 0')
                    elif letra == "c":
                        u_a = np.zeros(N)
                        plt.plot(eixo_t,u_a,color=(0.8,0,0.6,0.2), label ='T = 0')

                else:
                    if letra == "a":
                        if e<5:
                            alpha = e*0.1+0.1
                            plt.plot(eixo_t,Us[Fepkins-1],color = (alpha,0,1,alpha),label='T = {}' .format(e/10))
                        elif e<10:
                            alpha = e*0.1
                            plt.plot(eixo_t,Us[Fepkins-1],color = (0.4,alpha,1-alpha,alpha),label='T = {}' .format(e/10))
                    elif letra == "b":
                        alpha = 0.1*e
                        plt.plot(eixo_t,Us[Fepkins-1],color = (1,0,1,alpha),label='T = {}' .format(e/10))

                    elif letra == "c":
                        alpha = 0.1*e
                        plt.plot(eixo_t,Us[Fepkins-1],color = (0.8,0,0.3,alpha),label='T = {}' .format(e/10))

            if letra == "a":
                        plt.plot(eixo_t,Sol_final,color=(0,1,0,1),label ='T = 1')
            elif letra == "b":
                        plt.plot(eixo_t,Sol_final,color=(0.8,0,0.6,1), label ='T = 1')
            elif letra == "c":
                        plt.plot(eixo_t,Sol_final,color=(0.8,0,0.3,1), label ='T = 1')
            if metodo == "e":
                Metodo = "Euler"
            elif metodo == "c":
                Metodo = "Crank-Nicolson"
            plt.ylabel('Temperatura')
            plt.xlabel('Ponto X da barra')
            plt.title('Gráfico da Tarefa {}, letra {})\nAproximação pelo método de {}\n N = {}, λ = {}' .format(prime,letra,Metodo,N-1,Lambda))
            plt.xlim(0,1)
            plt.legend()
            plt.show()


        if letra != "c":
            Er = erro(USS,Us)
            print()
            print("----------------------------------")
            print("O maior erro foi de ", np.amax(Er))
            print()

        end = time.time()
        print("O número de passos foi ", Passos)
        print()
        print("o tempo total de processamento foi" , end - processa,"segundos")
        print()

        print()

        print("Digite 1 para reiniciar o programa ( 1 )")
        print("Digite 2 para sair ( 2 )")
        corona = int(input())


        if corona == 2:
            run_program = 0
            print("tchau")
main()

