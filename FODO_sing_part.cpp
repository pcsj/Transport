
#include "FODO_sing_part.h"


int dsMap(double L, double lunghezzatotale, int n_step)
{
	double ds = lunghezzatotale/n_step;
	int n_p = (int) floor(L/ds);
	return n_p;
}


#ifdef DEBUG
void prod(double * A, vector< vector <double> > N, double S)
#else
void prod(double * A, vector< vector <double> > N)
#endif
{
	double *H = new double[4];
	for (int i=0; i<4; i++) H[i] = 0.0;

	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			H[i] += A[j] * N[i][j];

	for (int i=0; i<4; i++)
		A[i]=H[i];

#ifdef DEBUG
	printf("S = %f\nx = %f\np_x = %f\ny = %f\np_y = %f\n",S,A[0],A[1],A[2],A[3]);
#endif
}


double *optics(vector< vector <double> > N, int i,bool * effettuato_con_successo)
{
	double omega, s, *A=new double[2];
	for (int j = 0; j < 2; j++) A[j] =0.;

	omega = acos((N[i][i]+N[i+1][i+1])*0.5);

	//omega=sqrt(-(N[i][i]-N[i+1][i+1])*(N[i][i]-N[i+1][i+1])+4*N[i+1][i]*N[i][i+1])/(N[i][i]+N[i+1][i+1]);
	s=sin(omega);
	
	if ( s*(N[i][i+1]) < 0.) s=-s;
	
	A[0]= (N[i][i] - cos(omega)) / s;
	A[1]= N[i][i+1] / s;

#ifdef DEBUG
	printf("\nN[%i][%i] = %f\nN[%i][%i] = %f",i,i,N[i][i],i+1,i+1,N[i+1][i+1]);
	printf("\nomega = %f",omega);
	printf("\nsin(omega) = %f",s);
	printf("\nA[0] = %f\nA[1] = %f",A[0],A[1]);
	printf("\n");
#endif

	*effettuato_con_successo=true;
	return A;
}


vector< vector <double> > prodo(vector< vector <double> > A, vector< vector <double> > B, int dimensione) 
{
	vector <vector <double> > H(4,vector<double>(4,0));

#ifdef DEBUG
	printf("\n");
	printf("prima ciclo prodo\n");
		for(int k=0; k < 4; k++)
		{
			printf("\n");
			for(int j=0; j < 4; j++) printf("%+12.8f  ", A[k][j]);
		}
	printf("\n");
		for(int k=0; k < 4; k++)
		{
			printf("\n");
			for(int j=0; j < 4; j++) printf("%+12.8f  ", B[k][j]);
		}
	printf("\n");
#endif

	for (int i=0; i<dimensione; i++)
		for (int j=0; j<dimensione; j++) 
        	for (int k=0; k<dimensione; k++)             
	        	        H[i][j] += A[i][k] * B[k][j];

	for (int a=0;a<dimensione;a++)
		for (int i=0;i<dimensione;i++)
			B[a][i]=H[a][i];
	
#ifdef DEBUG
	printf("\n");
	printf("dopo ciclo prodo\n");
		for(int k=0; k < 4; k++)
		{
			printf("\n");
			for(int j=0; j < 4; j++) printf("%+12.8f  ", B[k][j]);
		}
	printf("\n");
#endif

	return B;
}

#ifdef DEBUG
vector< vector <double> >  defocusing (vector< vector <double> > M, double d1, double S, FILE * file, int contatore)
#else
vector< vector <double> >  defocusing (vector< vector <double> > M, double d1, double S)
#endif
{
	M[0][0]=M[1][1]=cosh(d1*S);
	M[0][1]=sinh(d1*S)/d1;
	M[1][0]=sinh(d1*S)*d1;
	M[2][2]=M[3][3]=cos(d1*S);
	M[2][3]=sin(d1*S)/d1;
	M[3][2]=-sin(d1*S)*d1;

#ifdef DEBUG
	fprintf(file,"\nMATRICE DEFOC. dentro defocusing N # %d", contatore);
	for(int k=0; k < 4; k++)
	{
		fprintf(file,"\n");
		for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
	}
#endif

	return M;
}

#ifdef DEBUG
vector< vector <double> >  focusing (vector< vector <double> > M, double d1, double S, FILE * file, int contatore)
#else
vector< vector <double> >  focusing (vector< vector <double> > M, double d1, double S)
#endif
{
	M[0][0]=M[1][1]=cos(d1*S);
	M[0][1]=sin(d1*S)/d1;
	M[1][0]=-sin(d1*S)*d1;
	M[2][2]=M[3][3]=cosh(d1*S);
	M[2][3]=sinh(d1*S)/d1;
	M[3][2]=sinh(d1*S)*d1;	

#ifdef DEBUG
	fprintf(file,"\nMATRICE FOC. dentro focusing N # %d", contatore);		
	for(int k=0; k < 4; k++)
	{
		fprintf(file,"\n");
		for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
	}
#endif

	return M;
}


#ifdef DEBUG
vector< vector <double> >  drift (vector< vector <double> > M, double S, FILE * file, int contatore)
#else
vector< vector <double> >  drift (vector< vector <double> > M, double S)
#endif
{
	M[0][0]=M[1][1]=M[2][2]=M[3][3]=1.;
	M[0][1]=M[2][3]=S;

#ifdef DEBUG
	fprintf(file,"\nMATRICE DRIFT dentro drift N # %d", contatore);
	for(int k=0; k < 4; k++)
	{
		fprintf(file,"\n");
		for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
	}
#endif

	return M;
}


void scrivimatr2D (vector< vector <double> > M, FILE * output2)
{
	for(int k=0; k < 4; k++)
	{
		fprintf(output2,"\n");
		for(int j=0; j < 4; j++) fprintf(output2,"%+12.8f  ", M[k][j]);
	}
}


void scrividati (double s, double *A, double *B, FILE * file)
{
	fprintf(file,"\n %+7.4f",s);
	fprintf(file," %+10.5f",A[0]);			// alpha_x
	fprintf(file," %+10.5f",A[1]);			// beta_x
	fprintf(file," %+10.5f",B[0]);			// alpha_y
	fprintf(file," %+10.5f",B[1]);			// beta_y
}


void scrividati_ellissi (double s,double *AMaxMin, double *BMaxMin, FILE * file)
{
	fprintf(file,"\n %+7.4f",s);
	fprintf(file," %+10.5f",AMaxMin[1]);	// proiezione ellisse su asse x
	fprintf(file," %+10.5f",BMaxMin[1]);	// proiezione ellisse su asse p_x
	fprintf(file," %+10.5f",AMaxMin[0]);	// proiezione ellisse su asse y
	fprintf(file," %+10.5f",BMaxMin[0]);	// proiezione ellisse su asse p_y
}


void inizializza3D(double ***M, int a, int dimensione)
{
	M= new double**[a];
	for (int kk=0; kk < dimensione; kk++)
	{
		M[kk] = new double*[dimensione];
		for (int i=0; i < dimensione; i++) M[kk][i]= new double[dimensione];
	}
	
	for (int k = 0; k < a; k++)
		for (int i = 0; i < dimensione; i++)
			for (int j = 0; j<dimensione; j++)
				M[k][i][j]=0.0;
}


void inizializza2D(double **M, int dimensione)
{
	M= new double*[dimensione];
	for (int kk=0; kk < dimensione; kk++)
		M[kk] = new double[dimensione];
	
	for (int i = 0; i < dimensione; i++)
		for (int j = 0; j<dimensione; j++)
			M[i][j]=0.0;
}


void scrivi_pos_part(FILE * posizionePart,double *vett_i, double S)
{
	fprintf(posizionePart," %+10.5f",S);
	fprintf(posizionePart," %+10.5f",vett_i[0]);
	fprintf(posizionePart," %+10.5f",vett_i[1]);
	fprintf(posizionePart," %+10.5f",vett_i[2]);
	fprintf(posizionePart," %+10.5f",vett_i[3]);
	fprintf(posizionePart," %+10.5f",-(vett_i[0]));
	fprintf(posizionePart," %+10.5f\n",-(vett_i[2]));
}


vector< vector <double> > simil(vector< vector <double> > F,vector< vector <double> > OI, vector< vector <double> > O)
{
	vector <vector <double> > K(4,vector<double>(4,0.0));

	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			for (int a=0;a<4;a++)
				for (int l=0;l<4;l++) K[i][j]+=O[i][a]*F[a][l]*OI[l][j];

	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++) 
			F[i][j]=K[i][j];

	return F;
}


double * proiezioni_assi_ellissi(double *ottiche_x, double emittance)
{
	double * minimi_massimi = new double[2];
	for (int i = 0; i < 2; i++) minimi_massimi[i] = 0.0;
	minimi_massimi[0] = sqrt((ottiche_x[1])*emittance);
	minimi_massimi[1] = (ottiche_x[0] * emittance / (minimi_massimi[0]));
	return minimi_massimi;
}


void create_gnuplot_file_pos(string gnuplot_filename, string data_filename, double *lunghezza, int contatore, double estremo, double zmin, double zmax, string *keys)
{
	ofstream gnuplot_file;
	double lunghezza_percorsa=0.;
	gnuplot_file.open(gnuplot_filename.c_str());
	gnuplot_file << "#!/gnuplot" << endl;
	gnuplot_file << "FILE=\"" << data_filename << "\"" << endl;
	gnuplot_file << "set grid" << endl;
#if defined (CREATE_PNG)
	gnuplot_file << "set terminal png enhanced 15" << endl;
	gnuplot_file << "set output \"graph_" << keys[0] << ".png\"" << endl;
#elif defined (CREATE_EPS)
	gnuplot_file << "set terminal postscript eps enhanced colour solid rounded \"Helvetica\" 25" << endl;
	gnuplot_file << "set output \"graph_" << run_name << ".eps\"" << endl;
#endif
	gnuplot_file << "set xrange[" << zmin << ":" << zmax << "]" << endl;
	gnuplot_file << "set yrange[" << -estremo << ":" << estremo << "]" << endl;
	gnuplot_file << "set title  \""<< keys[1] <<" \"" << endl;
	gnuplot_file << "set xlabel \" " << keys[2] << "\"" << endl;
	gnuplot_file << "set ylabel \" " << keys[3] << "\"" << endl;
	for (int i = 0; i < contatore; i++)
	{
		lunghezza_percorsa+=lunghezza[i];
		gnuplot_file << "set arrow from " << lunghezza_percorsa<< "," << -estremo << " to "<< lunghezza_percorsa << ","<< estremo << " nohead lc rgb \"black\" lw 1" << endl;
	}
	gnuplot_file << "plot FILE u 1:2 w lines lt 1 lc rgb \"red\" lw 1 t \" " << keys[4] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:4 w lines lt 1 lc rgb \"blue\" lw 1 t \" " << keys[6] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:6 w lines lt 1 lc rgb \"red\" lw 1 ,\\" << endl;
	gnuplot_file << "FILE u 1:7 w lines lt 1 lc rgb \"blue\" lw 1 " << endl;

	//gnuplot_file << "FILE u 1:4 w lines lt 1 lc rgb \"orange\" lw 1 t \" " << keys[6] << "\",\\" << endl;
	//gnuplot_file << "FILE u 1:5 w lines lt 1 lc rgb \"dark-green\" lw 1 t \" " << keys[7] << "\"" << endl;
	gnuplot_file.close();
}

void create_gnuplot_file_opt(string gnuplot_filename, string data_filename, double *lunghezza, int contatore, double estremo, double zmin, double zmax, string *keys)
{
	ofstream gnuplot_file;
	double lunghezza_percorsa=0.;
	gnuplot_file.open(gnuplot_filename.c_str());
	gnuplot_file << "#!/gnuplot" << endl;
	gnuplot_file << "FILE=\"" << data_filename << "\"" << endl;
	gnuplot_file << "set grid" << endl;
#if defined (CREATE_PNG)
	gnuplot_file << "set terminal png enhanced 15" << endl;
	gnuplot_file << "set output \"graph_" << keys[0] << ".png\"" << endl;
#elif defined (CREATE_EPS)
	gnuplot_file << "set terminal postscript eps enhanced colour solid rounded \"Helvetica\" 25" << endl;
	gnuplot_file << "set output \"graph_" << run_name << ".eps\"" << endl;
#endif
	gnuplot_file << "set xrange[" << zmin << ":" << zmax << "]" << endl;
	gnuplot_file << "set yrange[" << -estremo << ":" << estremo << "]" << endl;
	gnuplot_file << "set title  \""<< keys[1] <<" \"" << endl;
	gnuplot_file << "set xlabel \" " << keys[2] << "\"" << endl;
	gnuplot_file << "set ylabel \" " << keys[3] << "\"" << endl;
	for (int i = 0; i < contatore; i++)
	{
		lunghezza_percorsa+=lunghezza[i];
		gnuplot_file << "set arrow from " << lunghezza_percorsa<< "," << -estremo << " to "<< lunghezza_percorsa << ","<< estremo << " nohead lc rgb \"black\" lw 1" << endl;
	}
	gnuplot_file << "plot FILE u 1:2 w lines lt 1 lc rgb \"red\" lw 1 t \" " << keys[4] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:3 w lines lt 1 lc rgb \"blue\" lw 1 t \" " << keys[5] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:4 w lines lt 1 lc rgb \"orange\" lw 1 t \" " << keys[6] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:5 w lines lt 1 lc rgb \"dark-green\" lw 1 t \" " << keys[7] << "\"" << endl;
	gnuplot_file.close();
}


#ifdef TEST_OPTICAL_FUNCTIONS
double *optics_T (double * A, int i, vector< vector <double> > O)
{
	double ottiche_x=0.;
	double ottiche_y=0.;
	ottiche_x = A[0];
	ottiche_y = A[1];

	A[0] = ottiche_x - ((O[i][i]) * (O[i+1][i]) * ottiche_y)  +  (2. * (O[i+1][i]) * (O[i][i+1]) * ottiche_x)  - ( (1./ottiche_y) * (O[i][i+1]) * (O[i+1][i+1]) * (1. + (ottiche_x * ottiche_x)));
	A[1]= ((O[i][i]) * (O[i][i]) * ottiche_y ) - (2. * (O[i][i]) * (O[i][i+1]) * ottiche_x)  + ( (1./ottiche_y) * (O[i][i+1]) * (O[i][i+1]) * (1. + (ottiche_x * ottiche_x)));
	return A;
}
#endif

void massimo_opt(double * optics_x , double * optics_y,double * massimo_temp)
{
	double max= *massimo_temp;
	max>fabs(optics_x[0])?max:max=fabs(optics_x[0]);
	max>fabs(optics_x[1])?max:max=fabs(optics_x[1]);
	max>fabs(optics_y[0])?max:max=fabs(optics_y[0]);
	max>fabs(optics_y[1])?max:max=fabs(optics_y[1]);
	* massimo_temp=max;
}


void massimo_pos(double * vett_i, double * massimo_temp)
{
	double max= *massimo_temp;
	max>fabs(vett_i[0])?max:max=fabs(vett_i[0]);
	max>fabs(vett_i[1])?max:max=fabs(vett_i[1]);
	max>fabs(vett_i[2])?max:max=fabs(vett_i[2]);
	max>fabs(vett_i[3])?max:max=fabs(vett_i[3]);
	* massimo_temp=max;
}


int main(int argc, char *argv[]) 
{ 
#ifdef __linux
	feenableexcept(2);
	feenableexcept(3);
#endif

	bool fallita_lettura_parametri = true;
	bool fallita_lettura_inputdistr = true;
	bool do_transport = false;
	bool do_optics = false;
	bool posso_fare_funzioni_ottiche = false;
	ifstream parametri;
	ifstream inputdistr;
	int nstep = 1;

	double gnuplot_ymax_opt=0.;
	bool calcola_ymax_opt = true;
#ifdef TEST_OPTICAL_FUNCTIONS
	double gnuplot_ymax_opt_T=0.;
	bool calcola_ymax_opt_T = true;
#endif
	double gnuplot_ymax_pos=0.;
	bool calcola_ymax_pos = true;
	double gnuplot_xmax_opt=0.;
	double gnuplot_xmax_pos=0.;
	bool calcola_ymax_ell = true;
	double gnuplot_ymax_ell=0.;


	double max_x=0.;
	double max_y=0.;
	bool coincidenza_xf_yf=false;
	bool confronto_pos=false;
	bool fail=false;
	double minimo_temp_x=100;
	double minimo_temp_y=100;
	double z_minimo_x=0.0;
	double z_minimo_y=1.0;

	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]) == "-p")
		{
			parametri.open(argv[i+1]);
			fallita_lettura_parametri=parametri.fail();
			i++;
		}
		else if (string(argv[i]) == "-i")
		{
			inputdistr.open(argv[i+1]);
			fallita_lettura_inputdistr=inputdistr.fail();
			i++;
		}
		else if (string(argv[i]) == "-transport")
		{
			do_transport=true;
		}
		else if (string(argv[i]) == "-xmax_opt")
		{
			gnuplot_xmax_opt=atof(argv[i+1]);
			i++;
		}
		else if (string(argv[i]) == "-xmax_pos")
		{
			gnuplot_xmax_pos=atof(argv[i+1]);
			i++;
		}
		else if (string(argv[i]) == "-ymax_opt")
		{
			gnuplot_ymax_opt=atof(argv[i+1]);
			calcola_ymax_opt = false;
			i++;
		}
#ifdef TEST_OPTICAL_FUNCTIONS
		else if (string(argv[i]) == "-ymax_opt_T")
		{
			gnuplot_ymax_opt_T=atof(argv[i+1]);
			calcola_ymax_opt_T = false;
			i++;
		}
#endif
		else if (string(argv[i]) == "-ymax_pos")
		{
			gnuplot_ymax_pos=atof(argv[i+1]);
			calcola_ymax_pos = false;
			i++;
		}
		else if (string(argv[i]) == "-ymax_ell")
		{
			gnuplot_ymax_ell=atof(argv[i+1]);
			calcola_ymax_ell = false;
			i++;
		}
		else if (string(argv[i]) == "-optimization")
		{
			confronto_pos=true;
		}
		else if (string(argv[i]) == "-optics")
		{
			do_optics=true;
		}
		else if (string(argv[i]) == "-nstep")
		{
			nstep = atoi(argv[i+1]);
			i++;
		}
		else
		{
			printf("Impossibile riconoscere il parametro %s\n",argv[i]);
		}
	}

	FILE * posizionePart=fopen("Posizione_Particelle.txt","w");

	FILE * ellissi=fopen("Parametri_Ellissi_Funz_Ottiche.txt","w");
	FILE * funzioni_ottiche=fopen("Funzioni_Ottiche.txt","w");
#ifdef TEST_OPTICAL_FUNCTIONS
	FILE * funzioni_ottiche_t=fopen("Funzioni_Ottiche_T.txt","w");
	FILE * ellissi_t=fopen("Parametri_Ellissi_Funz_Ottiche_T.txt","w");
#endif

#ifdef DEBUG
	FILE * matrici_iniziali=fopen("Matrici_Iniziali.txt","w");
	FILE * outputDEBUG=fopen("DEBUG.txt","w");
#endif

	string utile_per_contare;
	int conta_righe_parametri = 0;
	if (fallita_lettura_parametri || fallita_lettura_inputdistr)
	{
		printf("Impossibile aprire (o non definito) il file contenente i parametri\no il file contenente la distribuzione/particella iniziale\n");
		exit(204);
	}

// Cancelliamo i file che non ci servono
	if (!(do_transport))
	{
#if defined (__linux)
		system ("rm Posizione_Particelle.txt"); 
#elif defined (_WIN32) || defined (_WIN64)
		system ("del Posizione_Particelle.txt"); 
#endif		
	}



	double * dati_iniziali = new double[8];	// emittanza, energia, x, y, px, py
	for (int i = 0; i < 8; i++)
	{
		if(inputdistr.eof())
		{
			printf("Mancano dei dati iniziali!\n");
			exit(123);
		}
		inputdistr >> dati_iniziali[i];
	}
	inputdistr.clear();
	inputdistr.seekg(0,std::ios::beg);

	double emittanza = dati_iniziali[0];
	double energia = dati_iniziali[1];
	double *vett_i=new double[4];
	vett_i[0]=dati_iniziali[4];
	vett_i[1]=dati_iniziali[5];
	vett_i[2]=dati_iniziali[6];
	vett_i[3]=dati_iniziali[7];

	do
	{
		parametri >> utile_per_contare;
		if(parametri.eof()) break;
		parametri.ignore(1000, '\n');
		conta_righe_parametri++;
	}
	while(!parametri.eof());
	parametri.clear();
	parametri.seekg(0,std::ios::beg);

	// qui di leggono tutti i dati
	string *elemento=new string[conta_righe_parametri];
	double * lunghezza= new double[conta_righe_parametri];
	double * gradiente= new double[conta_righe_parametri];
	int contatore=0;
	for (int i = 0; i < conta_righe_parametri; i++)
	{
		parametri >> elemento[i];
		parametri >> gradiente[i];
		parametri >> lunghezza[i];
#ifdef DEBUG
		cout << "Tipo elemento: " << elemento[i] << ", gradiente: " << gradiente[i] << ", lunghezza: " << lunghezza[i] << endl;
#endif
		contatore++;
	}

#ifdef DEBUG
	printf("contatore: %d",contatore);
#endif

	vector <vector <double> > I(4,vector<double>(4,0.0));
	vector <vector <double> > K(4,vector<double>(4,0.0));
	vector <vector <double> > F(4,vector<double>(4,0.0));
	vector <vector <vector <double> > > Fx(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > Dx(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > OI(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > FxI(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > DxI(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > O(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));

	double *ottiche_x = new double[2];
	double *ottiche_y = new double[2];
	double *assi_ottiche_x = new double[2];
	double *assi_ottiche_y = new double[2];
	bool alpha_calcolato_con_successo=false;
	bool beta_calcolato_con_successo=false;

#ifdef TEST_OPTICAL_FUNCTIONS
	double *ottiche_x_t = new double[2];
	double *ottiche_y_t = new double[2];
	double *assi_ottiche_x_t = new double[2];
	double *assi_ottiche_y_t = new double[2];
	for (int i = 0; i < 2; i++) ottiche_x_t[i] = ottiche_y_t[i] = assi_ottiche_x_t[i] = assi_ottiche_y_t[i] = 0.;
#endif


	double gamma_beta=sqrt(2.0*energia/MP_MEV);
#ifdef cgs
	double gamma_v=gamma_beta*SPEED_OF_LIGHT*SPEED_OF_LIGHT*1.0e-4;
#else
	double gamma_v=gamma_beta*SPEED_OF_LIGHT;
#endif

	double *f1 =new double [contatore];
	double *d1 =new double [contatore];
	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=="F")
		{
			f1[i]=sqrt(gradiente[i]*CHARGE/(MP_KG*gamma_v));
#ifdef DEBUG
			cout<<"io sono f= "<<f1[i]*f1[i]<<endl;
#endif
		}
		if (elemento[i]=="D")
		{
			d1[i]=sqrt(gradiente[i]*CHARGE/(MP_KG*gamma_v));
#ifdef DEBUG
			cout<<"io sono d= "<<d1[i]*d1[i]<<endl;
#endif
		}
	}

#ifdef DEBUG
	for (int i=0;i<contatore;i++)
	{
		fprintf(outputDEBUG,"\ngrad. foc.    %+20.10f ", f1[i]*f1[i]);
		fprintf(outputDEBUG,"\ngrad. defoc.  %+20.10f ", d1[i]*d1[i]);
	}
#endif

	for (int i=0;i<contatore;i++)
	{
#ifdef DEBUG
		if (elemento[i]=="F")
			Fx[i]=focusing(Fx[i],f1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=="D")
			Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=="O")
			O[i]=drift(O[i],lunghezza[i],matrici_iniziali,i);
		else
			fprintf(outputDEBUG,"Elemento[%d] non riconosciuto\n", i);
#else
		if (elemento[i]=="F")
			Fx[i]=focusing(Fx[i],f1[i],lunghezza[i]);
		else if (elemento[i]=="D")
			Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i]);
		else if (elemento[i]=="O")
			O[i]=drift(O[i],lunghezza[i]);
#endif
	}

#ifdef DEBUG

	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=="O")
		{
//			if (O[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE DRIFT");
			scrivimatr2D(O[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}
		else if (elemento[i]=="F")
		{
//			if (Fx[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE FOC.");
			scrivimatr2D(Fx[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}
		else if (elemento[i]=="D")
		{
//			if (Dx[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE DEFOC.");
			scrivimatr2D(Dx[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}

	}

#endif	

/************************************************************************/
	
	if (do_optics)
	{
		vector <vector <double> > compos(4,vector<double>(4,0));

		if (elemento[0]=="O")
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=O[0][k][j];
		}
		else if (elemento[0]=="F")
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=Fx[0][k][j];
		}
		else if (elemento[0]=="D")
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=Dx[0][k][j];
		}

		for (int i=1;i<contatore;i++)
		{
			if (elemento[i]=="O")
				compos=prodo(O[i],compos,4);
			else if (elemento[i]=="F")
				compos=prodo(Fx[i],compos,4);
			else if (elemento[i]=="D")
				compos=prodo(Dx[i],compos,4);
		}
	
		for (int i=0;i<4;i++)
			for(int a=0;a<4;a++)
				F[i][a]=compos[i][a];

//		Calcolo Funzioni OTTICHE

		for (int i = 0; i < 2; i++) ottiche_x[i] = ottiche_y[i] = assi_ottiche_x[i] = assi_ottiche_y[i] = 0.;

		if ( (fabs((F[FOC][FOC]+F[FOC+1][FOC+1])*0.5) <= 1.) && (fabs((F[DEFOC][DEFOC]+F[DEFOC+1][DEFOC+1])*0.5) <= 1.))
			posso_fare_funzioni_ottiche = true;
		else cout << "Impossibile calcolare le funzioni ottiche!" << endl;

		if (posso_fare_funzioni_ottiche)
		{
			ottiche_x=optics(F,FOC,&alpha_calcolato_con_successo);
			ottiche_y=optics(F,DEFOC,&beta_calcolato_con_successo);
			assi_ottiche_x = proiezioni_assi_ellissi(ottiche_x, emittanza);
			assi_ottiche_y = proiezioni_assi_ellissi(ottiche_y, emittanza);

			if (calcola_ymax_opt) massimo_opt(ottiche_x,ottiche_y,&gnuplot_ymax_opt);
			if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
			if (calcola_ymax_ell) massimo_opt(assi_ottiche_x,assi_ottiche_y,&gnuplot_ymax_ell);
			if (alpha_calcolato_con_successo&&beta_calcolato_con_successo)
			{
				fprintf(funzioni_ottiche,"\n#%7c",'S');
				fprintf(funzioni_ottiche,"%10.8s","Alpha x");
				fprintf(funzioni_ottiche,"%10.7s","Beta x");
				fprintf(funzioni_ottiche,"%12.8s","Alpha y");
				fprintf(funzioni_ottiche,"%10.7s","Beta y");
				fprintf(ellissi,"%10s","x");
				fprintf(ellissi,"%11s","p_x");
				fprintf(ellissi,"%11s","y");
				fprintf(ellissi,"%11s","p_y");
			}

			scrividati(0.0,ottiche_x,ottiche_y,funzioni_ottiche);
			scrividati_ellissi(0.0,assi_ottiche_x,assi_ottiche_y,ellissi);


#ifdef TEST_OPTICAL_FUNCTIONS
			for (int i=0;i<2;i++)
			{
				if (fai_da_te_x&&fai_da_te_y)
				{
					ottiche_x_t[i]=paramIniz_X[i];
					ottiche_y_t[i]=paramIniz_Y[i];
				}
				else if (fai_da_te_x)
				{
					ottiche_x_t[i]=paramIniz_X[i];
					ottiche_y_t[i]=paramIniz_X[i];

				}
				else if (fai_da_te_y)
				{
					ottiche_x_t[i]=paramIniz_Y[i];
					ottiche_y_t[i]=paramIniz_Y[i];
				}
				else
				{
					ottiche_x_t[i]=ottiche_x[i];
					ottiche_y_t[i]=ottiche_y[i];
				}
			}
			assi_ottiche_x_t = proiezioni_assi_ellissi(ottiche_x_t, emittanza);
			assi_ottiche_y_t = proiezioni_assi_ellissi(ottiche_y_t, emittanza);
			scrividati(0.0,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
			scrividati_ellissi(0.0,assi_ottiche_x_t,assi_ottiche_y_t,ellissi_t);
			if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
#endif
		}
	}

	if(do_transport)
	{
		fprintf(posizionePart," %+10.5f",0.0);
		fprintf(posizionePart," %+10.5f",vett_i[0]);
		fprintf(posizionePart," %+10.5f",vett_i[1]);
		fprintf(posizionePart," %+10.5f",vett_i[2]);
		fprintf(posizionePart," %+10.5f",vett_i[3]);
		fprintf(posizionePart," %+10.5f",-vett_i[0]);
		fprintf(posizionePart," %+10.5f\n",-vett_i[2]);
	}
#ifdef DEBUG
		fprintf(outputDEBUG, "\nFODO:");
		scrivimatr2D(F,outputDEBUG);
#endif

/************************************************************************/

//ora primi dell'iterazione mi calcolo le micromappe Li di lunghezza S=L/n

	double lunghezzatotale=0.;
	for (int i = 0 ;i < contatore;i++)
		lunghezzatotale+=lunghezza[i];

#ifdef DEBUG
	for (int i=0; i < contatore; i++)
	{
		fprintf(outputDEBUG,"\n#step in elemento %d = %d",i, dsMap(lunghezza[i],lunghezzatotale,nstep));
	}
	fprintf(outputDEBUG,"\n");
#endif

	double S = 0.;
	
//Calcolo MICROMAPPE per il Drift/Focus/Defoc

	for (int i=0; i < contatore; i++)
	{
		S=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,nstep);
#ifdef DEBUG
		if (elemento[i] == "O")
		{
			O[i]=drift(O[i],S,matrici_iniziali,i);
			OI[i]=drift(OI[i],-S,matrici_iniziali,i);
		}
		else if (elemento[i] == "F")
		{
			Fx[i]=focusing(Fx[i],f1[i],S,matrici_iniziali,i);
			FxI[i]=focusing(FxI[i],f1[i],-S,matrici_iniziali,i);
		}
		else if (elemento[i] == "D")
		{
			Dx[i]=defocusing(Dx[i],d1[i],S,matrici_iniziali,i);
			DxI[i]=defocusing(DxI[i],d1[i],-S,matrici_iniziali,i);
		}
#else
		if (elemento[i] == "O")
		{
			O[i]=drift(O[i],S);
			OI[i]=drift(OI[i],-S);
		}
		else if (elemento[i] == "F")
		{
			Fx[i]=focusing(Fx[i],f1[i],S);
			FxI[i]=focusing(FxI[i],f1[i],-S);
		}
		else if (elemento[i] == "D")
		{
			Dx[i]=defocusing(Dx[i],d1[i],S);
			DxI[i]=defocusing(DxI[i],d1[i],-S);
		}
#endif
	}

/***********************************************************************/

double vett_i_temp_x,vett_i_temp_y;
double dl=0.;
double lunghezza_accumulata=0.0;

if(confronto_pos)
{
	for (int i=0;i<contatore;i++)
	{
		dl=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,nstep);
		if (elemento[i]=="O")
		{
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			{
				vett_i_temp_x=vett_i[0];
				vett_i_temp_y=vett_i[2];
#ifdef DEBUG
				prod(vett_i,O[i],S);
#else
				prod(vett_i,O[i]);
#endif
				if((vett_i[0]*vett_i_temp_x)<0.)
					z_minimo_x=S;

				if((vett_i[2]*vett_i_temp_y)<0.)
					z_minimo_y=S;
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
		if (elemento[i]=="F")
		{
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			{
				vett_i_temp_x=vett_i[0];
				vett_i_temp_y=vett_i[2];
#ifdef DEBUG
				prod(vett_i,Fx[i],S);
#else
				prod(vett_i,Fx[i]);
#endif
				if((vett_i[0]*vett_i_temp_x)<0.)
					z_minimo_x=S;

				if((vett_i[2]*vett_i_temp_y)<0.)
					z_minimo_y=S;
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}	
		if (elemento[i]=="D")
		{
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			{
				vett_i_temp_x=vett_i[0];
				vett_i_temp_y=vett_i[2];

#ifdef DEBUG				
				prod(vett_i,Dx[i],S);
#else
				prod(vett_i,Dx[i]);
#endif
				if((vett_i[0]*vett_i_temp_x)<0.)
					z_minimo_x=S;

				if((vett_i[2]*vett_i_temp_y)<0.)
					z_minimo_y=S;
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
	}

#ifdef DEBUG
	cout << "io sono z_min_x= "<<z_minimo_x<< "   io sono z_min_y= "<<z_minimo_y<<endl;
	cout << "io sono la differenza = "<< (fabs(z_minimo_x-z_minimo_y))<<endl;
#endif

	if (fabs(z_minimo_x-z_minimo_y)<(SHIFT))
		coincidenza_xf_yf=true;
	if((coincidenza_xf_yf==false)&&(confronto_pos==true))
		fail=true;
	
	dl=0.;
	lunghezza_accumulata=0.0;
	S=0.;
	for (int i=0;i<4;i++)
	{
		vett_i[i]=dati_iniziali[i+4];
	}
}

	for (int i=0;i<contatore;i++)
	{
		dl=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,nstep);
		if (elemento[i]=="O")
		{
#ifdef DEBUG
			fprintf(matrici_iniziali,"\n#Drift #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Drift #%d, dl = %f",i,dl);
#endif
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			{
#ifdef DEBUG
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
#endif
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,O[i],S);
#else
					prod(vett_i,O[i]);
#endif
					if ((confronto_pos==true)&&((vett_i[2]>0.01)&&(vett_i[0]>0.01)))
					{
						fail=true;
						continue;
					}

					if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (do_optics)
				{
					F=simil(F,OI[i],O[i]);
					if (posso_fare_funzioni_ottiche)
					{	
					ottiche_x=optics(F,FOC,&alpha_calcolato_con_successo);
					ottiche_y=optics(F,DEFOC,&beta_calcolato_con_successo);
					assi_ottiche_x = proiezioni_assi_ellissi(ottiche_x, emittanza);
					assi_ottiche_y = proiezioni_assi_ellissi(ottiche_y, emittanza);
					scrividati(S,ottiche_x,ottiche_y,funzioni_ottiche);
					scrividati_ellissi(S,assi_ottiche_x,assi_ottiche_y,ellissi);
					if (calcola_ymax_ell) massimo_opt(assi_ottiche_x,assi_ottiche_y,&gnuplot_ymax_ell);
					if (calcola_ymax_opt) massimo_opt(ottiche_x,ottiche_y,&gnuplot_ymax_opt);
#ifdef TEST_OPTICAL_FUNCTIONS
					ottiche_x_t=optics_T(ottiche_x_t,FOC,O[i]);
					ottiche_y_t=optics_T(ottiche_y_t,DEFOC,O[i]);
					assi_ottiche_x_t = proiezioni_assi_ellissi(ottiche_x_t, emittanza);
					assi_ottiche_y_t = proiezioni_assi_ellissi(ottiche_y_t, emittanza);
					scrividati(S,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
					scrividati_ellissi(S,assi_ottiche_x_t,assi_ottiche_y_t,ellissi_t);	
					if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
#endif
					}
				}
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
		else if (elemento[i]=="F")
		{
#ifdef DEBUG
			fprintf(matrici_iniziali,"\n#Drift #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Drift #%d, dl = %f",i,dl);
#endif
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			{
#ifdef DEBUG
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
#endif
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,Fx[i],S);
#else
					prod(vett_i,Fx[i]);
#endif
					if ((confronto_pos==true)&&((vett_i[2]>0.01)&&(vett_i[0]>0.01)))
					{
						fail=true;
						continue;
					}

					if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (do_optics)
				{
					if (posso_fare_funzioni_ottiche)
					{
						F=simil(F,FxI[i],Fx[i]);
						ottiche_x=optics(F,FOC,&alpha_calcolato_con_successo);
						ottiche_y=optics(F,DEFOC,&beta_calcolato_con_successo);		
						assi_ottiche_x = proiezioni_assi_ellissi(ottiche_x, emittanza);
						assi_ottiche_y = proiezioni_assi_ellissi(ottiche_y, emittanza);
						scrividati(S,ottiche_x,ottiche_y,funzioni_ottiche);
						scrividati_ellissi(S,assi_ottiche_x,assi_ottiche_y,ellissi);
						if (calcola_ymax_ell) massimo_opt(assi_ottiche_x,assi_ottiche_y,&gnuplot_ymax_ell);
						if (calcola_ymax_opt) massimo_opt(ottiche_x,ottiche_y,&gnuplot_ymax_opt);
#ifdef TEST_OPTICAL_FUNCTIONS
						ottiche_x_t=optics_T(ottiche_x_t,FOC,Fx[i]);
						ottiche_y_t=optics_T(ottiche_y_t,DEFOC,Fx[i]);
						assi_ottiche_x_t = proiezioni_assi_ellissi(ottiche_x_t, emittanza);
						assi_ottiche_y_t = proiezioni_assi_ellissi(ottiche_y_t, emittanza);
						scrividati(S,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
						scrividati_ellissi(S,assi_ottiche_x_t,assi_ottiche_y_t,ellissi_t);
						if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
#endif
					}
				}
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
		else if (elemento[i]=="D")
		{
#ifdef DEBUG
			fprintf(matrici_iniziali,"\n#Drift #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Drift #%d, dl = %f",i,dl);
#endif
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			{
#ifdef DEBUG
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
#endif
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,Dx[i],S);
#else
					prod(vett_i,Dx[i]);
#endif
					if ((confronto_pos==true)&&((vett_i[2]>0.01)&&(vett_i[0]>0.01)))
					{
						fail=true;
						continue;
					}

					if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (do_optics)
				{
					if (posso_fare_funzioni_ottiche)
					{
						F=simil(F,DxI[i],Dx[i]);
						ottiche_x=optics(F,FOC,&alpha_calcolato_con_successo);
						ottiche_y=optics(F,DEFOC,&beta_calcolato_con_successo);
						assi_ottiche_x = proiezioni_assi_ellissi(ottiche_x, emittanza);
						assi_ottiche_y = proiezioni_assi_ellissi(ottiche_y, emittanza);
						scrividati(S,ottiche_x,ottiche_y,funzioni_ottiche);
						scrividati_ellissi(S,assi_ottiche_x,assi_ottiche_y,ellissi);
						if (calcola_ymax_ell) massimo_opt(assi_ottiche_x,assi_ottiche_y,&gnuplot_ymax_ell);
						if (calcola_ymax_opt) massimo_opt(ottiche_x,ottiche_y,&gnuplot_ymax_opt);
#ifdef TEST_OPTICAL_FUNCTIONS
						ottiche_x_t=optics_T(ottiche_x_t,FOC,Dx[i]);
						ottiche_y_t=optics_T(ottiche_y_t,DEFOC,Dx[i]);
						assi_ottiche_x_t = proiezioni_assi_ellissi(ottiche_x_t, emittanza);
						assi_ottiche_y_t = proiezioni_assi_ellissi(ottiche_y_t, emittanza);
						scrividati(S,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
						scrividati_ellissi(S,assi_ottiche_x_t,assi_ottiche_y_t,ellissi_t);
						if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
#endif
						}
					}
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
	}

	fclose(funzioni_ottiche);
	fclose(posizionePart);
	fclose(ellissi);
	parametri.close();
	inputdistr.close();

#ifdef TEST_OPTICAL_FUNCTIONS
	fclose(funzioni_ottiche_t);
	fclose(ellissi_t);
#endif

	string *etichette_posizione = new string[8];
	string *etichette_ottiche = new string[8];
	string *etichette_ellissi = new string[8];
#ifdef TEST_OPTICAL_FUNCTIONS
	string *etichette_ottiche_T = new string[8];
	string *etichette_ellissi_T = new string[8];
#endif

	etichette_posizione[0] = "Posizione_Particelle";
	etichette_posizione[1] = "Posizione Particelle";
#ifdef cgs
	etichette_posizione[2] = "z (cm)";
	etichette_posizione[3] = "x/y (cm), p_x/p_y (g cm/s)";
#else 
	etichette_posizione[2] = "z (m)";
	etichette_posizione[3] = "x/y (m), p_x/p_y (Kg m/s)";
#endif
	etichette_posizione[4] = "x";
	etichette_posizione[5] = "p_x";
	etichette_posizione[6] = "y";
	etichette_posizione[7] = "p_y";

	etichette_ottiche[0] = "Funzioni_Ottiche";
	etichette_ottiche[1] = "Funzioni Ottiche";
#ifdef cgs
	etichette_ottiche[2] = "z (cm)";
#else
	etichette_ottiche[2] = "z (m)";
#endif
#if defined (CREATE_EPS)
	etichette_ottiche[3] = "{/Symbol a}, {/Symbol b}";
	etichette_ottiche[4] = "{/Symbol a}_x";
	etichette_ottiche[5] = "{/Symbol a}_y";
	etichette_ottiche[6] = "{/Symbol b}_x";
	etichette_ottiche[7] = "{/Symbol b}_y";
#else
#ifdef cgs
	etichette_ottiche[3] = "Alpha, Beta (cm)";
#else
	etichette_ottiche[3] = "Alpha, Beta (m)";
#endif
	etichette_ottiche[4] = "Alpha_x";
	etichette_ottiche[5] = "Beta_x";
	etichette_ottiche[6] = "Alpha_y";
	etichette_ottiche[7] = "Beta_y";
#endif

	etichette_ellissi[0] = "Parametri_Ellissi_Funz_Ottiche";
	etichette_ellissi[1] = "Parametri Ellissi Funz Ottiche";
#ifdef cgs
	etichette_ellissi[2] = "z (cm)";
	etichette_ellissi[3] = "X (cm), P (g cm/s)";
#else
	etichette_ellissi[2] = "z (m)";
	etichette_ellissi[3] = "X (m), P (Kg m/s)";
#endif
	etichette_ellissi[4] = "Xmax";
	etichette_ellissi[5] = "Pmax_x";
	etichette_ellissi[6] = "Ymax";
	etichette_ellissi[7] = "Pmax_y";

#ifdef TEST_OPTICAL_FUNCTIONS
	etichette_ottiche_T[0] = "Funzioni_Ottiche_T";
	etichette_ottiche_T[1] = "Funzioni ottiche test";
#ifdef cgs
	etichette_ottiche_T[3] = "z (cm)";
#else
	etichette_ottiche_T[3] = "z (m)";
#endif

#ifdef CREATE_EPS
	etichette_ottiche_T[3] = "{/Symbol a}, {/Symbol b}";
	etichette_ottiche_T[4] = "{/Symbol a}_x";
	etichette_ottiche_T[5] = "{/Symbol a}_y";
	etichette_ottiche_T[6] = "{/Symbol b}_x";
	etichette_ottiche_T[7] = "{/Symbol b}_y";
#else
#ifdef cgs
	etichette_ottiche_T[3] = "Alpha, Beta (cm)";
#else
	etichette_ottiche_T[3] = "Alpha, Beta (m)";
#endif
	etichette_ottiche_T[4] = "Alpha_x";
	etichette_ottiche_T[5] = "Beta_x";
	etichette_ottiche_T[6] = "Alpha_y";
	etichette_ottiche_T[7] = "Beta_y";
#endif

	etichette_ellissi_T[0] = "Parametri_Ellissi_Funz_Ottiche_T";
	etichette_ellissi_T[1] = "Parametri Ellissi Funz Ottiche_T";
#ifdef cgs
	etichette_ellissi_T[2] = "z (cm)";
	etichette_ellissi_T[3] = "X (cm), P (g cm/s)";
#else
	etichette_ellissi_T[2] = "z (m)";
	etichette_ellissi_T[3] = "X (m), P (Kg m/s)";
#endif
	etichette_ellissi_T[4] = "Xmax";
	etichette_ellissi_T[5] = "Pmax_x";
	etichette_ellissi_T[6] = "Ymax";
	etichette_ellissi_T[7] = "Pmax_y";
#endif

/***********************************************************/

if (fail==false)
{
	if (do_transport)
	{
		if (gnuplot_xmax_pos > 0.) 
		{
			create_gnuplot_file_pos( "Posizione.plt", "Posizione_Particelle.txt", lunghezza, contatore, gnuplot_ymax_pos ,0.0, gnuplot_xmax_pos, etichette_posizione);
			//create_gnuplot_file( "Posizione_thick.plt", "Posizione_Particelle_thick.txt", lunghezza, contatore, gnuplot_ymax_pos ,0.0, gnuplot_xmax_pos, etichette_posizione);
		}
		else 
		{
			create_gnuplot_file_pos( "Posizione.plt", "Posizione_Particelle.txt", lunghezza, contatore, gnuplot_ymax_pos ,0.0, lunghezza_accumulata, etichette_posizione);
			//create_gnuplot_file( "Posizione_thick.plt", "Posizione_Particelle_thick.txt", lunghezza, contatore, gnuplot_ymax_pos ,0.0, lunghezza_accumulata, etichette_posizione);		
		}
		system ("gnuplot Posizione.plt");
		//system ("gnuplot Posizione_thick.plt");
	}
	else exit(203);

	if ((alpha_calcolato_con_successo&&beta_calcolato_con_successo)&&(do_optics))
	{
		if (gnuplot_xmax_opt > 0.) create_gnuplot_file_opt( "Funzioni_Ottiche.plt", "Funzioni_Ottiche.txt", lunghezza, contatore, gnuplot_ymax_opt ,0.0, gnuplot_xmax_opt, etichette_ottiche);
		else create_gnuplot_file_opt( "Funzioni_Ottiche.plt", "Funzioni_Ottiche.txt", lunghezza, contatore, gnuplot_ymax_opt ,0.0, lunghezza_accumulata, etichette_ottiche);
		create_gnuplot_file_opt( "Parametri_Ellissi.plt", "Parametri_Ellissi_Funz_Ottiche.txt", lunghezza, contatore, gnuplot_ymax_ell ,0.0, lunghezza_accumulata, etichette_ellissi);
		system ("gnuplot Funzioni_Ottiche.plt");
		system ("gnuplot Parametri_Ellissi.plt");

#ifdef TEST_OPTICAL_FUNCTIONS
		create_gnuplot_file_opt( "Funzioni_Ottiche_T.plt", "Funzioni_Ottiche_T.txt", lunghezza, contatore, gnuplot_ymax_opt ,0.0, lunghezza_accumulata, etichette_ottiche_T);
		create_gnuplot_file_opt( "Parametri_Ellissi_T.plt", "Parametri_Ellissi_Funz_Ottiche_T.txt", lunghezza, contatore, gnuplot_ymax_ell ,0.0, lunghezza_accumulata, etichette_ellissi_T);
		system ("gnuplot Funzioni_Ottiche_T.plt");
		system ("gnuplot Parametri_Ellissi_T.plt");
#endif
	}
	else
	{
#if defined (__linux)
		system ("rm Funzioni_Ottiche.txt"); 
		system ("rm Parametri_Ellissi_Funz_Ottiche.txt");
#ifdef TEST_OPTICAL_FUNCTIONS
		system ("rm Funzioni_Ottiche_T.txt");
		system ("rm Parametri_Ellissi_Funz_Ottiche_T.txt");
#endif
#elif defined (_WIN32) || defined (_WIN64)
		system ("del Funzioni_Ottiche.txt"); 
		system ("del Parametri_Ellissi_Funz_Ottiche.txt"); 
#ifdef TEST_OPTICAL_FUNCTIONS
		system ("del Funzioni_Ottiche_T.txt");
		system ("del Parametri_Ellissi_Funz_Ottiche.txt");
#endif
#endif
	} 
	}
else
{
#if defined (__linux)
		system ("rm Posizione_Particelle.txt");
		system ("rm Funzioni_Ottiche.txt"); 
		system ("rm Parametri_Ellissi_Funz_Ottiche.txt");
#ifdef TEST_OPTICAL_FUNCTIONS
		system ("rm Funzioni_Ottiche_T.txt");
		system ("rm Parametri_Ellissi_Funz_Ottiche_T.txt");
#endif
#elif defined (_WIN32) || defined (_WIN64)
		system ("del Posizione_Particelle.txt"); 
		system ("del Funzioni_Ottiche.txt"); 
		system ("del Parametri_Ellissi_Funz_Ottiche.txt");
#ifdef TEST_OPTICAL_FUNCTIONS
		system ("del Funzioni_Ottiche_T.txt");
		system ("del Parametri_Ellissi_Funz_Ottiche.txt");
#endif
#endif
		exit(204);
}	

#ifdef DEBUG
	fclose(outputDEBUG);
	fclose(matrici_iniziali);
#endif

}

