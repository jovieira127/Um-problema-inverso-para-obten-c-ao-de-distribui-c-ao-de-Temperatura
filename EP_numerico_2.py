
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

'''Função que cria a matriz'''
def crie_matriz(n_linhas, n_colunas, valor):
    matriz = []
    for i in range(n_linhas):      
        linha = []
        for j in range(n_colunas):
            linha.append(valor)
        matriz.append(linha)
    return matriz

'''Funções para o cálculo da temperatura'''
def funcao_r(t):
    f = 10*(1+math.cos(5*t))
    return f

def gh(h):
    return (1/h)

def funcaoC(x,t,p,h):
    if x >= p-(h/2) and x <= p+h/2:
        f = gh(h)*funcao_r(t)
    else:
        f = 0
    return f 

'''fatorar e resolver matriz'''
def fatora_Matriz(A, B):
    U = []
    U.append(A[0])
    L = []
    n = len (A)
    
    for i in range(n-1):
        l = B[i]/U[i]
        L.append(l)  
        d = A[i+1]-L[i]**2*U[i]
        U.append(d)
            
    return L, U

def Resolve(A, B,lista):
    L, U = fatora_Matriz(A,B)
    n = len(U)
    z=[]
    y=[]
    x=[]
    z.append(lista[0])
    for k in range(1, len(A)):
        z1 = lista[k] - L[k-1]*z[k-1]
        z.append(z1)
    for j in range (len(z)):
        y1 = z[j]/U[j]
        y.append(y1)
    x.append(y[(len(y)-1)])
    for i in range (len(A)-2, -1, -1):
        x1 = y[i]-L[i]*x[0]
        x.insert(0, x1)
    return x

''' Tarefa A método de crank'''
def crank(N,M,lambida,p,h,intervalot):

    Resultados = crie_matriz(int(M+1),int(N+1),0.0)
    '''criando diagonal'''
    diagonal =[]
    for i in range (int(N)-1):
        diagonal.append(1+lambida)
    '''criando subdiagonal'''
    subdiagonal =[]
    for i in range (int(N)-1):
        subdiagonal.append(-lambida/2)
    for k in range (int(M)):
        b=[]
        for i in range (1, int(N)):
            bi = Resultados[k][i] + (lambida/2)*(Resultados[k][i-1]- 2*Resultados[k][i]+Resultados[k][i+1]) +(intervalot/2)*(funcaoC(i*intervalox,k*intervalot,p,h)+funcaoC(i*intervalox,(k+1)*intervalot,p,h))
            b.append(bi)
        c=Resolve(diagonal, subdiagonal,b)
        for j in range (1, int(N)):
            Resultados[k+1][j] = c[j-1]

    return Resultados

'''Copiar a ultima linha da matriz para para outra matriz'''
def copiaultimalinha(matriz1,N,matrizfinal,posicao):
    for i in range (N+1):
        matrizfinal[posicao][i] = matriz1[N][i]
    return matrizfinal    

'''Tarefa B Calcular produtos internos e montar a matriz de produtos internos'''
def produtointernovetormatriz(matriz,vetor,N,nf,resultado):

    for i in range(nf):
        produto = 0
        for z in range (N+1):
            produto += matriz[i][z]*vetor[z]
        resultado[i] = produto    

    return resultado

def produtointernomatrizmatriz(matriz,produtointerno,N,nf):

    for i in range(nf):
        for z in range (nf):
            produto = 0
            for j in range(N+1):
                 produto += matriz[i][j]*matriz[z][j]
            produtointerno[i][z] = produto    

    return produtointerno

'''construcao do gráfico'''
def grafico(matriz, intervalo_x, intervalo_t):
    interv = 0.1/intervalo_t
    cont = 0
    while cont < len(matriz):
        posicao = 0
        x = []
        y = []
        for j in range(len(matriz[0])):
            x.append(posicao)
            y.append(matriz[int(cont)][j])
            posicao += intervalo_x  
        plt.plot(x, y, label = cont/(interv*10))
        plt.legend()
        plt.grid()
        plt.xlabel("Posição")
        plt.ylabel("Temperatura")
        plt.title("Temperatura x Distância com N = %i e Lambida = %.3f" %((1/intervalo_x), ((intervalo_t/(intervalo_x**2)))))
        cont += interv      
    return plt.show()

'''funcoes auxiliar que foi utilizada para verificar se o produto interno estava correto'''
def imprimematriz(Matriz, Tamanho):
    for i in range(Tamanho):
        for j in range(Tamanho):
            print(int(Matriz[i][j]),end=" ")
        print("\n")
        
def imprimevetor(Vetor, Tamanho):
    for i in range(Tamanho):
        print(Vetor[i])
        
        
        
'''Tarefa C, construindo o método LDLt e resolvendo a matriz'''
def fatoraLU(A):  
    U = np.copy(A)  
    n = np.shape(U)[0]  
    L = np.eye(n)  
    for j in np.arange(n-1):  
        for i in np.arange(j+1,n):  
            L[i,j] = U[i,j]/U[j,j]  
            for k in np.arange(j+1,n):  
                U[i,k] = U[i,k] - L[i,j]*U[j,k]  
            U[i,j] = 0
    return L, U

def Resolve_LU(A,b):
    L, U = fatoraLU(A)
    n = len(U[0])
    z=[]
    x=[]
    z.append(b[0])
    for i in range(1, n):
        z1 = 0
        for j in range(0, i):
            z1 += L[i][j]*z[j]
        z.append((b[i] - z1)/(L[i][i]))
    
    x.append((z[n-1])/(U[n-1][n-1]))
    for i in range (n-2, -1, -1):
        xa = 0
        tamx = len(x)-1
        for j in range (n-1, i, -1):
            xa += U[i][j]*x[tamx]
            tamx-=1
        x.insert(0, (z[i] - xa)/(U[i][i]))    
    return x
    
CASO = int(input("Digite 1 para o teste a ou 2 para o teste b, nao fizemos os testes c e d\n"))
                    
T = 1
'''Teste A'''
if CASO == 1:
    '''parametros iniciais para calcular u1,u2, etc'''
    N = 128
    intervalox = 1/N
    M = N
    intervalot = float(T/M)
    lambida = float(N)
    Resultados = crie_matriz(int(M+1),int(N+1),0.0)

    '''para calcular os produtos internos'''
    p = []
    p.append(0.35)
    h = intervalox 
    nf = 1
    Prodint = crie_matriz(nf,nf,0.0)
    iguala = []
    iguala.append(0)
    matrizfinal = crie_matriz(int(nf),int(N+1),0.0)
    
    '''crank, aqui a tarefa A é executada'''
    for i in range(nf):
        Resultados = crank(N,M,lambida,p[i],h,intervalot)
        matrizfinal = copiaultimalinha(Resultados,N,matrizfinal,i)
        '''para provar que o método de crank esta funcionando'''
        '''grafico(Resultados, intervalox, intervalot)'''
    '''matrizfinal eh a matriz com todos os resultados obtidos por crank, um em cada linha'''
    
    '''Aproximacao eh o vetor UT'''
    Aproximacao = []
    for i in range (int(N+1)):
        Aproximacao.append(7*matrizfinal[0][i])

    '''Aqui a tarefa B é executada'''       
    Prodint = produtointernomatrizmatriz(matrizfinal,Prodint,N,nf)       
    iguala = produtointernovetormatriz(matrizfinal,Aproximacao,N,nf,iguala)

    '''Imprimindo a matriz de produtos internos, foram colocados apenas os inteiros para facilitar a visualizacao'''
    print("matriz de inteiros de produtos internos =")
    imprimematriz(Prodint, nf)
    
    '''Imprimindo o vetor de resultados de produtos internos'''
    print("vetor de produtos internos =")
    imprimevetor(iguala, nf)
    print("\n")

    '''Aqui a Tarefa C eh executada'''
    o = Resolve_LU(Prodint,iguala)
    print("coeficientes", o)
    
    
if CASO == 2:
    '''parametros iniciais para calcular u1,u2, etc'''
    N = 128
    intervalox = 1/N
    M = N
    intervalot = float(T/M)
    lambida = float(N)
    Resultados = crie_matriz(int(M+1),int(N+1),0.0)
    
    '''para calcular os produtos internos'''
    p = []
    p.append(0.15)
    p.append(0.3)
    p.append(0.7)
    p.append(0.8)
    h = intervalox 
    nf = 4
    Prodint = crie_matriz(nf,nf,0.0)
    posicao = 0
    matrizfinal = crie_matriz(int(nf),int(N+1),0.0)
    
    '''crank, aqui a tarefa A é executada'''
    for i in range(nf):
        Resultados = crank(N,M,lambida,p[i],h,intervalot)
        matrizfinal = copiaultimalinha(Resultados,N,matrizfinal,i)
        '''para provar que o método de crank esta funcionando'''
        '''grafico(Resultados, intervalox, intervalot)'''
    '''matrizfinal eh a matriz com todos os resultados obtidos por crank, um em cada linha'''
    
    iguala = []
    for i in range (nf):
        iguala.append(0)
        
    '''Aproximacao eh o vetor UT'''    
    Aproximacao = []
    for i in range (int(N+1)):
        Aproximacao.append(2.3*matrizfinal[0][i]+3.7*matrizfinal[1][i]+0.3*matrizfinal[2][i]+4.2*matrizfinal[3][i])

    '''matriz de respostas'''
    Prodint = produtointernomatrizmatriz(matrizfinal,Prodint,N,nf)       

    '''vetor de resultados'''
    iguala = produtointernovetormatriz(matrizfinal,Aproximacao,N,nf,iguala)

    
    '''Imprimindo a matriz de produtos internos, foram colocados apenas os inteiros para facilitar a visualizacao'''
    print("matriz de de inteiros produtos internos =")
    imprimematriz(Prodint, nf)
    
    '''Imprimindo o vetor de resultados de produtos internos'''
    print("vetor de produtos internos =")
    imprimevetor(iguala, nf)
    print("\n")

    
    o = Resolve_LU(Prodint,iguala)
    print("coeficientes", o)



