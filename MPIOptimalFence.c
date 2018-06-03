#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "ConvexHull.h"
#include <mpi.h>
//#include <omp.h>

#define DMaxArboles 	27
#define DMaximoCoste 999999
# define MTAG1 1
# define MTAG2 2
#define NUM_THREADS 16

  //////////////////////////
 // Estructuras de datos //
//////////////////////////


// Tree structure definition.
struct  Arbol
{
	int	  IdArbol;
	Point Coord;			// tree position
	int Valor;				// Value
	int Longitud;			// wood quantity
};
typedef struct Arbol TArbol, *PtrArbol;



// Forest structure definition.
struct Bosque
{
	int 		NumArboles;
	TArbol 	Arboles[DMaxArboles];
};
typedef struct Bosque TBosque, *PtrBosque;



// Problem respresentation
struct ListaArboles
{
	int 		NumArboles;
 	float		Coste;
	float		CosteArbolesCortados;
	float		CosteArbolesRestantes;
	float		LongitudCerca;
	float		MaderaSobrante;
	int 		Arboles[DMaxArboles];
};
typedef struct ListaArboles TListaArboles, *PtrListaArboles;


// static coordinates.
typedef Point TVectorCoordenadas[DMaxArboles], *PtrVectorCoordenadas;



typedef enum {false, true} bool;


  ////////////////////////
 // Variables Globales //
////////////////////////

TBosque ArbolesEntrada;



  //////////////////////////
 // Funtion definition   //
//////////////////////////
bool LeerFicheroEntrada(char *PathFicIn); //Input file to read
bool GenerarFicheroSalida(TListaArboles optimo, char *PathFicOut);//Generate the output file
bool CalcularCercaOptima(PtrListaArboles Optimo);//Calculate the optimal fence
void OrdenarArboles(); // Sort the Trees
bool CalcularCombinacionOptima(int PrimeraCombinacion, int UltimaCombinacion, PtrListaArboles Optimo);//Calculate the optimal combination
int EvaluarCombinacionListaArboles(int Combinacion); //Evaluate the combination of the list trees
int ConvertirCombinacionToArboles(int Combinacion, PtrListaArboles CombinacionArboles);//Converts the combination to trees
int ConvertirCombinacionToArbolesTalados(int Combinacion, PtrListaArboles CombinacionArbolesTalados);//convert the combination to trees choped
void ObtenerListaCoordenadasArboles(TListaArboles CombinacionArboles, TVectorCoordenadas Coordenadas);// Obtain the trees's cordinate lists 
float CalcularLongitudCerca(TVectorCoordenadas CoordenadasCerca, int SizeCerca);// Calculate the length of the fence
float CalcularDistancia(int x1, int y1, int x2, int y2);//Calculate the distance
int CalcularMaderaArbolesTalados(TListaArboles CombinacionArboles);//Calculate the wood of the trees already cut
int CalcularCosteCombinacion(TListaArboles CombinacionArboles);//calculate the cost of the combination
void MostrarArboles(TListaArboles CombinacionArboles);//Show the trees



int main(int argc, char *argv[])
{
	if (argc<2 || argc>3)
		printf("Error Argumentos");


	int myid , numprocs , i , islave ;
	
	MPI_Init (&argc ,&argv);
	MPI_Status status;
	MPI_Comm_size ( MPI_COMM_WORLD ,& numprocs );
	MPI_Comm_rank ( MPI_COMM_WORLD ,& myid );

	

	TListaArboles Optimo;
	char *posicion;
	int namesrclen;
	unsigned long len;
	char dest[255],outfile[255];
	double tpivot1=0,tpivot2=0; //time counting
	struct timeval tim;
	//omp_set_num_threads(4);
	
	if ( myid == 0) {
		gettimeofday(&tim, NULL);
    	tpivot1 = tim.tv_sec+(tim.tv_usec/1000000.0);
	}

	//Capture first token time
    if (!LeerFicheroEntrada(argv[1]))
	{
		printf("Error opening input file\n");
		exit(1);
	}

	if (!CalcularCercaOptima(&Optimo))
	{
		printf("Error computing the optimal solution\n");
		exit(1);
	}

	if ( myid == 0) {
		if (argc==2)
		{	
			len=strlen(argv[1]);
			posicion = strchr(argv[1], '.');
			memset(dest, '\0', sizeof(dest));
			namesrclen=posicion-argv[1];
			strncpy(dest,argv[1],namesrclen);
			sprintf(outfile,"%s.res", dest);

			if (!GenerarFicheroSalida(Optimo, outfile))
			{
				printf("Error writing output file\n");
				exit(1);
			}
		}
		else
		{
			if (!GenerarFicheroSalida(Optimo, argv[2]))
			{
				printf("Error writing output file\n");
				exit(1);
			}
		}
		
		gettimeofday(&tim, NULL);	
	    tpivot2 = (tim.tv_sec+(tim.tv_usec/1000000.0));
	    printf("\nTime %.6lf\n", tpivot2-tpivot1);
		//printf("\n Number of threads used : %d\n",omp_get_max_threads());
		//exit(0);
	}
	MPI_Finalize ();
}



bool LeerFicheroEntrada(char *PathFicIn)//Input file to read
{
	FILE *FicIn;
	int a;

	FicIn=fopen(PathFicIn,"r");
	if (FicIn==NULL)
	{
		perror("Opening input file");
		return false;
	}
	printf("Input data:\n");

	// Reading the number of trees in the input woods.
	if (fscanf(FicIn, "%d", &(ArbolesEntrada.NumArboles))<1)
	{
		perror("Reading input woods");
		return false;
	}
	printf("\tTrees: %d.\n",ArbolesEntrada.NumArboles);

	// Reading tree atributes.
	for(a=0;a<ArbolesEntrada.NumArboles;a++)
	{
		ArbolesEntrada.Arboles[a].IdArbol=a+1;
		// Reading  x, y, value, Longitud.
		if (fscanf(FicIn, "%d %d %d %d",&(ArbolesEntrada.Arboles[a].Coord.x), &(ArbolesEntrada.Arboles[a].Coord.y), &(ArbolesEntrada.Arboles[a].Valor), &(ArbolesEntrada.Arboles[a].Longitud))<4)
		{
			perror("Reading tree info");
			return false;
		}
		printf("\tTree %d-> (%d,%d) Value:%d, Long:%d.\n",a+1,ArbolesEntrada.Arboles[a].Coord.x, ArbolesEntrada.Arboles[a].Coord.y, ArbolesEntrada.Arboles[a].Valor, ArbolesEntrada.Arboles[a].Longitud);
	}
	
	return true;
}



bool GenerarFicheroSalida(TListaArboles Optimo, char *PathFicOut)//Generate the output file
{
	FILE *FicOut;
	int a;

	FicOut=fopen(PathFicOut,"w+");
	if (FicOut==NULL)
	{
		perror("Opening output file.");
		return false;
	}

	// Writting the number of trees to cut.
	if (fprintf(FicOut, "Cutting %3d trees           : ", Optimo.NumArboles)<1)
	{
		perror("writing number of cutting trees");
		return false;
	}

	for(a=0;a<Optimo.NumArboles;a++)
	{
		// Writing the id of a tree to cut.
		if (fprintf(FicOut, "%3d ",ArbolesEntrada.Arboles[Optimo.Arboles[a]].IdArbol)<1)
		{
			perror("writing id of a tree to cut.");
			return false;
		}
	}

	// Writing the total missing wood
	if (fprintf(FicOut, "\nFence long./missing Wood    : \t%6.2f (%6.2f)",   Optimo.LongitudCerca, Optimo.MaderaSobrante)<1)
	{
		perror("writing total missing wood");
		return false;
	}

	// Writing the cost of trees to cut
	if (fprintf(FicOut, "\nValue of the cutted Trees   : \t%6.2f", Optimo.CosteArbolesCortados)<1)
	{
		perror("writing value of the cutted Trees");
		return false;
	}

	// Writing the value of the remaining trees
	if (fprintf(FicOut, "\nValue of the Remaining Trees: \t%6.2f\n", Optimo.CosteArbolesRestantes)<1)
	{
		perror("writing value of the Remaining Trees");
		return false;
	}
	
	return true;
}


bool CalcularCercaOptima(PtrListaArboles Optimo)//Calculate the optimal fence
{
	int MaxCombinaciones;

	/* Computing total combinations */
	MaxCombinaciones = (int) pow(2.0,ArbolesEntrada.NumArboles);

	// Sort trees in increasing order by x,y
	OrdenarArboles();

	/* Computing optimal */
	Optimo->NumArboles = 0;
	Optimo->Coste = DMaximoCoste;
	CalcularCombinacionOptima(1, MaxCombinaciones, Optimo);

	return true;
}


// 
void OrdenarArboles()// Sort the Trees
{
  	
  	int a,b;
  	//TArbol aux;
	double start_time, run_time;
	start_time = MPI_Wtime(); 
	//#pragma omp parallel for default(none) shared(ArbolesEntrada) private(a,b,aux)
	for(a=0; a<(ArbolesEntrada.NumArboles-1); a++)
	{
		for(b=a; b<ArbolesEntrada.NumArboles; b++)
		{
			if ( ArbolesEntrada.Arboles[b].Coord.x < ArbolesEntrada.Arboles[a].Coord.x ||
				 (ArbolesEntrada.Arboles[b].Coord.x == ArbolesEntrada.Arboles[a].Coord.x && ArbolesEntrada.Arboles[b].Coord.y < ArbolesEntrada.Arboles[a].Coord.y) )
			{
				
				TArbol aux;
				// aux=a
				aux.Coord.x = ArbolesEntrada.Arboles[a].Coord.x;
				aux.Coord.y = ArbolesEntrada.Arboles[a].Coord.y;
				aux.IdArbol = ArbolesEntrada.Arboles[a].IdArbol;
				aux.Valor = ArbolesEntrada.Arboles[a].Valor;
				aux.Longitud = ArbolesEntrada.Arboles[a].Longitud;

				// a=b
				ArbolesEntrada.Arboles[a].Coord.x = ArbolesEntrada.Arboles[b].Coord.x;
				ArbolesEntrada.Arboles[a].Coord.y = ArbolesEntrada.Arboles[b].Coord.y;
				ArbolesEntrada.Arboles[a].IdArbol = ArbolesEntrada.Arboles[b].IdArbol;
				ArbolesEntrada.Arboles[a].Valor = ArbolesEntrada.Arboles[b].Valor;
				ArbolesEntrada.Arboles[a].Longitud = ArbolesEntrada.Arboles[b].Longitud;

				// b=aux
				ArbolesEntrada.Arboles[b].Coord.x = aux.Coord.x;
				ArbolesEntrada.Arboles[b].Coord.y = aux.Coord.y;
				ArbolesEntrada.Arboles[b].IdArbol = aux.IdArbol;
				ArbolesEntrada.Arboles[b].Valor = aux.Valor;
				ArbolesEntrada.Arboles[b].Longitud = aux.Longitud;
			}
		}
	}
	run_time = MPI_Wtime()- start_time;
	printf("run rime %f seconds \n",run_time);
}




// Computing the optimal combination in a range
bool CalcularCombinacionOptima(int PrimeraCombinacion, int UltimaCombinacion, PtrListaArboles Optimo)//Calculate the optimal combination
{
	int myid , numprocs , islave ;
	printf("Ultima combination %d\n", UltimaCombinacion );
	MPI_Status status;
	MPI_Comm_size ( MPI_COMM_WORLD ,&numprocs );
	MPI_Comm_rank ( MPI_COMM_WORLD ,&myid );

	int Combinacion, MejorCombinacion=0, CosteMejorCombinacion; //mejor=better
	int Coste;
	int Send[2], SendCoste, SendCombinacion;

	TListaArboles CombinacionArboles;
	TVectorCoordenadas CoordArboles, CercaArboles;
	int NumArboles, PuntosCerca;//puntos=point
	float MaderaArbolesTalados;

  	//printf("Evaluating combinations: \n");
	char buffer[100]; // buffer to store the data(MPI_Pack/MPI_Unpack) 
	int position; // the position of the buffer
	
	CosteMejorCombinacion = Optimo->Coste;
	double start_time, run_time;
        start_time = MPI_Wtime();

   	//MPI_Bcast(&UltimaCombinacion,1,MPI_INTEGER, 0,MPI_COMM_WORLD);
    //int nThreads=0;
    //omp_set_num_threads(NUM_THREADS);
	//#pragma omp parallel default(none) private(Combinacion, Coste) shared(nThreads, myid, PrimeraCombinacion, UltimaCombinacion, numprocs,CosteMejorCombinacion, MejorCombinacion, Send) 
//{
	//#pragma omp master
	//nThreads=omp_get_num_threads();
	
	//#pragma omp for
	for (Combinacion=myid+PrimeraCombinacion; Combinacion<UltimaCombinacion; Combinacion+=numprocs)
	{
    	
    	//printf("\tC%d -> \t",Combinacion);
		Coste = EvaluarCombinacionListaArboles(Combinacion);
		if ( Coste < CosteMejorCombinacion )
		{
			//#pragma omp critical
			//if (Coste < CosteMejorCombinacion)
			//{
				CosteMejorCombinacion = Coste;
				MejorCombinacion = Combinacion;
				Send[0] =CosteMejorCombinacion;
				Send[1] =MejorCombinacion;
			//}
			
		}
	}
//}
/*    if  (nThreads == NUM_THREADS) {  
      	printf("%d OpenMP threads were used.\n", NUM_THREADS);    
   	}  
   	else {  
      	printf("Expected %d OpenMP threads, but %d were used.\n", NUM_THREADS, nThreads);  
   	}*/

	printf("Coste of rank %d is %d\n",myid, CosteMejorCombinacion);
	//printf("Combinacion of rank %d is %d\n",myid, MejorCombinacion);
	
	if ( myid != 0) 
	{
		// Send the contents of buffer to processors 
		MPI_Send ((void*)Send, 2, MPI_INTEGER, 0 , MTAG2 , MPI_COMM_WORLD);
	} 
	else 
	{
		for ( islave =1; islave < numprocs ; islave ++) {
			MPI_Recv ((void*)Send,2 ,MPI_INTEGER, islave, MTAG2 , MPI_COMM_WORLD , &status );
					//Unpack the contents of buffer 
			printf("Receiving Coste %d from rank %d\n",Send[1], islave);
			SendCoste=Send[0];
			SendCombinacion=Send[1];

			if ( SendCoste < CosteMejorCombinacion )
			{
				CosteMejorCombinacion = SendCoste;
				MejorCombinacion = SendCombinacion;
	      	  //   printf("***");
			}
		}
		run_time = MPI_Wtime() - start_time;
        printf("run time CalcularCombinacionOptima %f seconds\n",run_time);

		if (CosteMejorCombinacion == Optimo->Coste)
			return false;  // No se ha encontrado una combinacin mejor. A better combination has not been found.

		// Asignar combinacin encontrada.Assign combination found.
		ConvertirCombinacionToArbolesTalados(MejorCombinacion, Optimo);
		Optimo->Coste = CosteMejorCombinacion;
		// Calcular estadisticas óptimo. Calculate optimal statistics.
		NumArboles = ConvertirCombinacionToArboles(MejorCombinacion, &CombinacionArboles);
		ObtenerListaCoordenadasArboles(CombinacionArboles, CoordArboles);
		PuntosCerca = chainHull_2D( CoordArboles, NumArboles, CercaArboles );

		Optimo->LongitudCerca = CalcularLongitudCerca(CercaArboles, PuntosCerca);
		MaderaArbolesTalados = CalcularMaderaArbolesTalados(*Optimo);
		Optimo->MaderaSobrante = MaderaArbolesTalados - Optimo->LongitudCerca;
		Optimo->CosteArbolesCortados = CosteMejorCombinacion;
		Optimo->CosteArbolesRestantes = CalcularCosteCombinacion(CombinacionArboles);

		return true;
	}

}



int EvaluarCombinacionListaArboles(int Combinacion)//Evaluate the combination of the list trees
{
	TVectorCoordenadas CoordArboles, CercaArboles;
	TListaArboles CombinacionArboles, CombinacionArbolesTalados;
	int NumArboles, NumArbolesTalados, PuntosCerca, CosteCombinacion;
	float LongitudCerca, MaderaArbolesTalados;

	// Convertimos la combinacin al vector de arboles no talados.We convert the combination to the vector of uncut trees.
	NumArboles = ConvertirCombinacionToArboles(Combinacion, &CombinacionArboles);
  
	// Obtener el vector de coordenadas de arboles no talados.Get the coordinate vector of trees not cut down.
	ObtenerListaCoordenadasArboles(CombinacionArboles, CoordArboles);

	// Calcular la cerca Calculate the fence
	PuntosCerca = chainHull_2D( CoordArboles, NumArboles, CercaArboles );

	/* Evaluar si obtenemos suficientes �boles para construir la cerca Evaluate if we get enough trees to build the fence */
	LongitudCerca = CalcularLongitudCerca(CercaArboles, PuntosCerca);

	// Evaluar la madera obtenida mediante los arboles talados.Evaluate the wood obtained through the felled trees.
	// Convertimos la combinacin al vector de arboles no talados.We convert the combination to the vector of uncut trees.
	NumArbolesTalados = ConvertirCombinacionToArbolesTalados(Combinacion, &CombinacionArbolesTalados);
  //printf(" %d arboles cortados: ",NumArbolesTalados);
  //MostrarArboles(CombinacionArbolesTalados);
  MaderaArbolesTalados = CalcularMaderaArbolesTalados(CombinacionArbolesTalados);
  //printf("  Madera:%4.2f  \tCerca:%4.2f ",MaderaArbolesTalados, LongitudCerca);
	if (LongitudCerca > MaderaArbolesTalados)
	{
		// Los arboles cortados no tienen suficiente madera para construir la cerca.The cut trees do not have enough wood to build the fence.
    //printf("\tCoste:%d",DMaximoCoste);
    return DMaximoCoste;
	}

	// Evaluar el coste de los arboles talados. Evaluate the cost of felled trees.
	CosteCombinacion = CalcularCosteCombinacion(CombinacionArbolesTalados);
  //printf("\tCoste:%d",CosteCombinacion);
  
	return CosteCombinacion;
}


int ConvertirCombinacionToArboles(int Combinacion, PtrListaArboles CombinacionArboles)//Converts the combination to trees
{
	int arbol=0;

	CombinacionArboles->NumArboles=0;
	CombinacionArboles->Coste=0;

	while (arbol<ArbolesEntrada.NumArboles)
	{
		if ((Combinacion%2)==0)
		{
			CombinacionArboles->Arboles[CombinacionArboles->NumArboles]=arbol;
			CombinacionArboles->NumArboles++;
			CombinacionArboles->Coste+= ArbolesEntrada.Arboles[arbol].Valor;
		}
		arbol++;
		Combinacion = Combinacion>>1;
	}

	return CombinacionArboles->NumArboles;
}


int ConvertirCombinacionToArbolesTalados(int Combinacion, PtrListaArboles CombinacionArbolesTalados)//convert the combination to trees choped
{
	int arbol=0;

	CombinacionArbolesTalados->NumArboles=0;
	CombinacionArbolesTalados->Coste=0;

	while (arbol<ArbolesEntrada.NumArboles)
	{
		if ((Combinacion%2)==1)
		{
			CombinacionArbolesTalados->Arboles[CombinacionArbolesTalados->NumArboles]=arbol;
			CombinacionArbolesTalados->NumArboles++;
			CombinacionArbolesTalados->Coste+= ArbolesEntrada.Arboles[arbol].Valor;
		}
		arbol++;
		Combinacion = Combinacion>>1;
	}

	return CombinacionArbolesTalados->NumArboles;
}



void ObtenerListaCoordenadasArboles(TListaArboles CombinacionArboles, TVectorCoordenadas Coordenadas)// Obtain the trees's cordinate lists 
{
	int c, arbol;

	for (c=0;c<CombinacionArboles.NumArboles;c++)
	{
    	arbol=CombinacionArboles.Arboles[c];
		Coordenadas[c].x = ArbolesEntrada.Arboles[arbol].Coord.x;
		Coordenadas[c].y = ArbolesEntrada.Arboles[arbol].Coord.y;
	}
}


	
float CalcularLongitudCerca(TVectorCoordenadas CoordenadasCerca, int SizeCerca)// Calculate the length of the fence
{
	int x;
	float coste;
	
	//#pragma omp parallel for reduction(+:coste)
	for (x=0;x<(SizeCerca-1);x++)
	{
		coste+= CalcularDistancia(CoordenadasCerca[x].x, CoordenadasCerca[x].y, CoordenadasCerca[x+1].x, CoordenadasCerca[x+1].y);
	}
	return coste;
}



float CalcularDistancia(int x1, int y1, int x2, int y2)//Calculate the distance
{
	return(sqrt(pow((double)abs(x2-x1),2.0)+pow((double)abs(y2-y1),2.0)));
}



int CalcularMaderaArbolesTalados(TListaArboles CombinacionArboles)//Calculate the wood of the trees already cut
{
	int a;
	int LongitudTotal=0;

	//#pragma omp parallel for reduction(+:LongitudTotal)
	for (a=0;a<CombinacionArboles.NumArboles;a++)
	{
		LongitudTotal += ArbolesEntrada.Arboles[CombinacionArboles.Arboles[a]].Longitud;
	}
	return(LongitudTotal);
}



int CalcularCosteCombinacion(TListaArboles CombinacionArboles)//calculate the cost of the combination
{
	int a;
	int CosteTotal=0;

	//#pragma omp parallel for reduction(+:CosteTotal)
	for (a=0;a<CombinacionArboles.NumArboles;a++)
	{
		CosteTotal += ArbolesEntrada.Arboles[CombinacionArboles.Arboles[a]].Valor;
	}
	return(CosteTotal);
}


void MostrarArboles(TListaArboles CombinacionArboles)//Show the trees
{
	int a;

	for (a=0;a<CombinacionArboles.NumArboles;a++)
		printf("%d ",ArbolesEntrada.Arboles[CombinacionArboles.Arboles[a]].IdArbol);

  	for (;a<ArbolesEntrada.NumArboles;a++)
    	printf("  ");  
}
