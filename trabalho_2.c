#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 100

typedef struct {
    double elementos[MAX][MAX];
    int ordem;
} matriz;

matriz LerMatriz();
matriz Submatriz(matriz, int);
double Cofator(matriz, int);
double Determinante(matriz);
void SistemaTriangularSuperior(matriz, double *, double *);
void SistemaTriangularInferior(matriz, double *, double *);
void DecomposicaoLU(matriz, double *, double *);

// Bibi e Rick

matriz LerMatriz() {
    matriz matriz;
    int i, j;

    do {
        printf("Defina a ordem da matriz (1-%d): ", MAX);
        scanf("%d", &matriz.ordem);
    }
    while (matriz.ordem < 1 || matriz.ordem > MAX);

    for (i = 0; i < matriz.ordem; i++)
    {
        for (j = 0; j < matriz.ordem; j++)
        {
            printf("Defina o elemento da %d%c linha e %d%c coluna: ", i + 1, 166, j + 1, 166);
            scanf("%lf", &matriz.elementos[i][j]);
        }
    }

    return matriz;
}

//Remove a última linha e uma coluna qualquer
//pelo índice (0 até matriz.ordem - 1)
matriz Submatriz(matriz a, int coluna) {
    int n = a.ordem - 1;
    int i, j;

    for (i = 0; i < n; i++)
    {
        for (j = coluna; j < n; j++)
        {
            a.elementos[i][j] = a.elementos[i][j + 1];
        }
    }

    a.ordem--;
    return a;
}

double Cofator(matriz a, int coluna) {
    return (1 + coluna + a.ordem) % 2 == 0
        ? Determinante(Submatriz(a, coluna))
        : -1 * Determinante(Submatriz(a, coluna));
}

//O i sempre é igual a matriz.ordem para ficar
//mais fácil calcular a submatriz
double Determinante(matriz a) {
  int j;
    if (a.ordem == 1) {
        return a.elementos[0][0];
    }
    else {
        double resultado = 0;
        int n = a.ordem - 1;

        for (j = 0; j < a.ordem; j++)
        {
            if (a.elementos[n][j] != 0) {
                resultado += a.elementos[n][j] * Cofator(a, j);
            }
        }

        return resultado;
    }
}

void SistemaTriangularSuperior(matriz a, double b[], double solucao[]) {
    int n = a.ordem - 1;
    int i, j;
    solucao[n] = b[n] / a.elementos[n][n];

    for (i = n - 1; i >= 0; i--)
    {
        double x = b[i];

        for (j = i + 1; j < a.ordem; j++)
        {
            x -= a.elementos[i][j] * solucao[j];
        }

        x /= a.elementos[i][i];
        solucao[i] = x;
    }
}

void SistemaTriangularInferior(matriz a, double b[], double solucao[]) {
  int i, j;
    solucao[0] = b[0] / a.elementos[0][0];

    for (i = 1; i < a.ordem; i++)
    {
        double x = b[i];

        for (j = 0; j < i; j++)
        {
            x -= a.elementos[i][j] * solucao[j];
        }

        x /= a.elementos[i][i];
        solucao[i] = x;
    }
}

void DecomposicaoLU(matriz a, double b[], double solucao[]) {
    matriz l, u;
    int i, h, j, k;
    l.ordem = u.ordem = a.ordem;

    for (j = 0; j < a.ordem; j++)
    {
        u.elementos[0][j] = a.elementos[0][j];
    }

    l.elementos[0][0] = 1;

    for (i = 1; i < a.ordem; i++)
    {
        l.elementos[i][i] = 1;
        l.elementos[i][0] = a.elementos[i][0] / u.elementos[0][0];
    }

    for (h = 1; h < a.ordem; h++)
    {
        for (j = h; j < a.ordem; j++)
        {
            u.elementos[h][j] = a.elementos[h][j];

            for (k = 0; k < h; k++)
            {
                u.elementos[h][j] -= l.elementos[h][k] * u.elementos[k][j];
            }

        }

        for (i = h + 1; i < a.ordem; i++)
        {
            l.elementos[i][h] = a.elementos[i][h];

            for (k = 0; k < h; k++)
            {
                l.elementos[i][h] -= l.elementos[i][k] * u.elementos[k][h];
            }

            l.elementos[i][h] /= u.elementos[h][h];
        }
    }

    double y[MAX];
    SistemaTriangularInferior(l, b, y);
    SistemaTriangularSuperior(u, y, solucao);
}

// Rick e Bibi
int simetrica(matriz a) {
  int i, j;
   for(i = 0; i < a.ordem; i++) {
      for(j = 0; j < a.ordem; j++) {
        if(i != j) {
          if(a.elementos[i][j] != a.elementos[j][i]){
            return 0;
          }
        }
      }
    }
    return 1;
}

int convergencia(matriz a) {
  int i;
  for(i = 1; i <= a.ordem; i++) {
     while(a.ordem > 0) {
       if(Determinante(a) != 0)
     		a.ordem--;
     	else
     		return 0;
    }
  }
  return 1;
}

int gaussJordan(matriz a, double b[], double x[]) {
  int i, j, k;
  double pivo, m;
  for(i = 0; i < a.ordem; i++)
    a.elementos[i][a.ordem] = b[i];

  for(k = 0; k < a.ordem; k++) {
    pivo = a.elementos[k][k];
    if(pivo == 0) return 0;
    for(i = 0; i < a.ordem; i++) {
      m = a.elementos[i][k] / pivo;
      for(j = 0; j < a.ordem + 1; j++) {
        if(i != k) {
          a.elementos[i][j] = a.elementos[i][j] - m * a.elementos[k][j];
        }
      }
    }
  }

  for(i = 0; i < a.ordem; i++) {
    for(j = 0; j < a.ordem + 1; j++) {
      printf("%.2lf          ", a.elementos[i][j]);
    }
    printf("\n");
  }

  for(i = 0; i < a.ordem; i++)
    for(j = 0; j < a.ordem + 1; j++) {
      if(i == j) {
        x[i] = a.elementos[i][a.ordem] / a.elementos[i][j];
      }
    }

  return 1;
}

int preCholesky(matriz a) {
  int flag = simetrica(a);
  int i;

  if(flag) {
    for(i = 1; i <= a.ordem; i++) {
       while(a.ordem > 0) {
    	if(Determinante(a) > 0)
    		a.ordem--;
    	else
    		return 0;
   		}
    }
    return 1;
  }
  return 0;
}

void transposta(matriz a, matriz *lT) {
  int i, j;
  for(i = 0; i < a.ordem; i++) {
    for(j = 0; j < a.ordem; j++) {
      (*lT).elementos[i][j] = a.elementos[j][i];
    }
  }
}

void Cholesky(matriz a, double b[], double solucao[]) {
  int i, j, k;
  double somatorio;
  matriz l, lT;
  l.ordem = lT.ordem = a.ordem;

  for(i = 0; i < a.ordem; i++) {
    for(j = 0; j < a.ordem; j++) {
      l.elementos[i][j] = 0;

      if(i == j) {
        somatorio = 0;
        for(k = 0; k <= i - 1; k++) {
          somatorio += pow(l.elementos[i][k], 2);
        }
        l.elementos[i][j] = sqrt(a.elementos[i][j] - somatorio);
      }
      else {
        somatorio = 0;

        for(k = 0; k <= j - 1; k++) {
          somatorio += l.elementos[i][k] * l.elementos[j][k];
        }

        if(j < i && l.elementos[j][j] != 0) {
        	l.elementos[i][j] = (a.elementos[i][j] - somatorio) / l.elementos[j][j];
		}
      }
    }
  }

  transposta(l, &lT);

  for(i = 0; i < a.ordem; i++) {
    for(j = 0; j < a.ordem; j++) {
      printf("%.2lf          ", l.elementos[i][j]);
    }
    printf("\n");
  }

  printf("\n\n");

  for(i = 0; i < a.ordem; i++) {
    for(j = 0; j < a.ordem; j++) {
      printf("%.2lf          ", lT.elementos[i][j]);
    }
    printf("\n");
  }

  double y[MAX];
  SistemaTriangularInferior(l, b, y);
  SistemaTriangularSuperior(lT, y, solucao);
}

int criterioLinhas(matriz a) {
	double s = 0;
	double maior = 0;
	int i, j;
	
	for(i = 0; i < a.ordem; i++) {
		s = 0;
		
		for(j = 0; j < a.ordem; j++) {
			if(i != j)
				s += fabs(a.elementos[i][j] / a.elementos[i][i]);
		}
		
		if(s > maior) {
			maior = s;	
		}
		
		if(maior >= 1) {
			return 0;
		}
	}
	return 1;
}

int criterioColunas(matriz a) {
	double s = 0;
	double maior = 0;
	int i, j;
	
	for(j = 0; j < a.ordem; j++) {
		s = 0;
		
		for(i = 0; i < a.ordem; i++) {
			if(i != j)
				s += fabs(a.elementos[i][j] / a.elementos[j][j]);
		}
		
		if(s > maior) {
			maior = s;	
		}
		
		if(maior >= 1) {
			return 0;
		}
	}
	return 1;
}

int criterioSassenfield(matriz a) {
	int i, j, n;
	double b[a.ordem], s1, s2;
	double maior = 0;
	for(i = 0; i < a.ordem; i++) {
		s1 = 0;
		s2 = 0;
		for(j = 0; j <= i - 1; j++) {
			s1 += fabs(a.elementos[i][j] / a.elementos[i][i]) * b[j];
		}
		for(j = i + 1; j < a.ordem; j++) {
			s2 += fabs(a.elementos[i][j] / a.elementos[i][i]);	
		}
		b[i] = s1 + s2;
		if(b[i] > maior) maior = b[i];
		if(maior >= 1) return 0;
	}
	return 1;
}

int preGaussSeidel(matriz a) {
  int i;
  
  for(i = 0; i < a.ordem; i++) {
  	if(a.elementos[i][i] == 0) return 0;
  }
  
  if(Determinante(a) == 0) return 0;
  
  if(!criterioSassenfield(a) && !criterioLinhas(a)) return 0;
  
  return 1;
}

int criterioJacobiESeidel(double x0[], double x[], double e, int ordem) {
	int i;
	double maiorNum = 0;
	double maiorDen = 0;
	double num, den;
	for(i = 0; i < ordem; i++) {
		num = fabs(x[i] - x0[i]);
		den = fabs(x[i]);
		if(num > maiorNum) maiorNum = num;
		if(den > maiorDen) maiorDen = den;
	}
	if(maiorDen != 0) return ((maiorNum / maiorDen) < e);
}

int gaussSeidel(matriz a, double b[], double e, double x0[], int max, double x[], int *it) {
	int i, j, k;
	double somatorio = 0;
	double aux[a.ordem];
	for(k = 0; k < max; k++) {
		(*it)++;
		for(i = 0; i < a.ordem; i++) {
			somatorio = 0;
			for(j = 0; j < a.ordem; j++) {
				if(i != j) {
					somatorio += a.elementos[i][j] * x0[j];
				}
			}
			x[i] = (b[i] - somatorio) / a.elementos[i][i];
			aux[i] = x0[i];
			x0[i] = x[i];
		}
		if(criterioJacobiESeidel(aux, x, e, a.ordem)) return 1;
	}
	return 0;
}

int preJacobi(matriz a) {
  int i;
  
  for(i = 0; i < a.ordem; i++) {
  	if(a.elementos[i][i] == 0) return 0;
  }
  
  if(Determinante(a) == 0) return 0;
  
  if(!criterioColunas(a) && !criterioLinhas(a)) return 0;
  
  return 1;
}

void mudaX0(double x[], double x0[], int ordem) {
	int i;
	for(i = 0; i < ordem; i++)
		x0[i] = x[i];
}

int Jacobi(matriz a, double b[], double e, double x0[], int max, double x[], int *it) {
	int i, j, k, l;
	double somatorio = 0;
	for(k = 0; k < max; k++) {
		(*it)++;
		for(i = 0; i < a.ordem; i++) {
			somatorio = 0;
			for(j = 0; j < a.ordem; j++) {
				if(i != j) {
					somatorio += a.elementos[i][j] * x0[j];
				}
			}
			x[i] = (b[i] - somatorio) / a.elementos[i][i];
		}
		if(criterioJacobiESeidel(x0, x, e, a.ordem)) return 1;
		
		mudaX0(x, x0, a.ordem);
	}
	return 0;
}

int main() {
    matriz A = LerMatriz();
    // matriz A;
    
//    A.ordem = 4;
//    A.elementos[0][0] = 4;
//    A.elementos[0][1] = -2;
//    A.elementos[0][2] = 1;
//    A.elementos[1][0] = 2;
//    A.elementos[1][1] = -6;
//    A.elementos[1][2] = -1;
//    A.elementos[2][0] = 1;
//    A.elementos[2][1] = 1;
//    A.elementos[2][2] = -3;

    double b[] = { -9, 20, 12, -4}, x[4];
    double e = 0.1;
    double x0[] = { 0, 0, 0, 0 };
    int max = 10;
    int it = 0;
    
    //for (int i = 0; i < matriz.ordem; i++)
    //{
    //    for (int j = 0; j < matriz.ordem; j++)
    //    {
    //        printf("%d ", matriz.elementos[i][j]);
    //    }

    //    printf("\n");
    //}

//    if(convergencia(A)) {
//      gaussJordan(A, b, x);
//      printf("%lf %lf %lf\n", x[0], x[1], x[2]);
//    }

//	if(preJacobi(A)) {
//		if(Jacobi(A, b, e, x0, max, x, &it)) {
//			printf("\n%.3lf %.3lf %.3lf %.3lf com %d iteracoes.", x[0], x[1], x[2], x[3], it); 
//		} else {
//			printf("\n O metodo nao pode ser resolvido com %d iteracoes.", it);
//		}
//	} else {
//		printf("\n O metodo nao converge.");
//	}

//	if(preGaussSeidel(A)) {
//		if(gaussSeidel(A, b, e, x0, max, x, &it)) {
//			printf("\n%lf %lf %lf %lf com %d iteracoes.", x[0], x[1], x[2], x[3], it); 
//		} else {
//			printf("\n O metodo nao pode ser resolvido com %d iteracoes.", it);
//		}
//	} else {
//		printf("\n O metodo nao converge.");
//	}
	
    // if(preCholesky(A)){
    //   Cholesky(A, b, x);
    //   printf("%lf %lf %lf %lf \n", x[0], x[1], x[2], x[3]);
    // }
    // DecomposicaoLU(matriz, b, x);
    // printf("%lf %lf %lf\n", x[0], x[1], x[2]);
}
