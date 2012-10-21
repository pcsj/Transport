$GRADIENTI_INIZIALI=50.0
$GRADIENTI_FINALI=150.0
$NUM_CICLI=10
$L_TOT=0.5
$L_IN=0.2
$L_ELEM=0.06
$ENERGY=30
$DRIFT_F=0.8

if (!(Test-Path "Plot")){
new-Item -type directory ("Plot")}

if (!(Test-Path "Input")){
new-Item -type directory ("Input")}

$numero_cicli=0

$grad_3=${GRADIENTI_INIZIALI}
while ($a -le $NUM_CICLI) #qui grad_3
{
    $grad_2=${GRADIENTI_INIZIALI}
    while ($b -le $NUM_CICLI) #qui grad_2
    {
        $grad_1=${GRADIENTI_INIZIALI}
        while ($c -le $NUM_CICLI) #qui grad_1
        {
            $L_1=0.05
            while ($d -le $NUM_CICLI) #qui L_1
            {
                $L_2=0.001
                while ($e -le $NUM_CICLI) #qui L_2
                {
                
                $L_3=$L_TOT-$L_1-$L_2
            
	            $riga="O 0.0 $L_IN"
	            out-file -filepath .\parametri.txt -inputobject $riga -encoding ASCII
    
                $riga="D $grad_1 $L_ELEM"
	            out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII

	            $riga="O 0.0 $L_2"
	            out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII

	            $riga="F $grad_1 $L_ELEM"
	            out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII
            
	            $riga="O 0.0 $L_1"
	            out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII
    
                $riga="F $grad_2 $L_ELEM"
	            out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII

	            $riga="O 0.0 $L_3"
	            out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII

	            $riga="D $grad_3 $DRIFT_F"
	            out-file -filepath .\parametri.txt -inputobject $riga -append -encoding ASCII            
    
                .\FODO_I.exe -p .\parametri.txt -i .\inputdata.txt -transport -nstep 200 -optimization -ymax_pos 0.3
    	    
                if (Test-Path ".\graph_Posizione_Particelle.png")
	            {
	               New-Item -Force -type directory ".\Plot\GRAD1_${grad_1}__GRADF2_${grad_2}__GRADD2_${grad_3}"
                   Move-Item -Force graph_Posizione_Particelle.png  .\Plot\GRAD1_${grad_1}__GRADF2_${grad_2}__GRADD2_${grad_3}.png
                   Move-Item -Force Posizione_Particelle.txt  .\Plot\GRAD1_${grad_1}__GRADF2_${grad_2}__GRADD2_${grad_3}.txt
                   Move-Item -Force Posizione.plt  .\Plot\Posizione_${numero_cicli}.plt                                              
                }
                $numero_cicli=1+$numero_cicli
                $L_2=$L_2+0.001
                $e=$e+1
            
            }
            $L_2=0.001
            $e=1    
            $L_1=$L_1+0.001
            $d=$d+1
            
            }
            $d=1
            $grad_1=$grad_1+1
            $c=$c+1
         }
         $c=1
         $grad_2=$grad_2+1
         $b=$b+1
    }
$b=1
$grad_3=$grad_3+1
$a=$a+1
}

            #if ($prova)
            #{
            #    while ($g -le $NUMBER_OF_STEPS_ENERGY)
            #    {
            #        $energy=$energy+2;
            #        echo $energy
	        #        $riga="1"
	        #        out-file -filepath .\inputdata.txt -inputobject $riga -encoding ASCII
	        #        $riga="$energy"
	        #        out-file -filepath .\inputdata.txt -inputobject $riga -append -encoding ASCII
	        #        $riga="0.0 0.0 0.05 0.05"
	        #        out-file -filepath .\inputdata.txt -inputobject $riga -append -encoding ASCII

            #        .\FODO_sing_part.exe -p .\parametri.txt -i .\inputdata.txt -optics -transport -nstep $NUMBER_OF_STEP_PER_SIM -compare_X $X_X_RIF $X_P_RIF -compare_Y $Y_X_RIF $Y_P_RIF -perc 0.1
                    
            #        if (Test-Path ".\graph_Posizione_Particelle.png")
	        #        {
            #           Move-Item -Force graph_Posizione_Particelle.png  .\Posizione_Particelle\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
            #        }
            #        if (Test-Path ".\graph_Funzioni_Ottiche.png")
            #        {            
	        #           Move-Item -Force graph_Funzioni_Ottiche.png  .\Funzioni_Ottiche\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
	        #        }
            #        if (Test-Path ".\graph_Funzioni_Ottiche_T.png")
            #        {
	        #           Move-Item -Force graph_Funzioni_Ottiche_T.png  .\Funzioni_Ottiche_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
	        #        }
            #        if (Test-Path ".\graph_Parametri_Ellissi_Funz_Ottiche.png")
	        #        {
            #           Move-Item -Force graph_Parametri_Ellissi_Funz_Ottiche.png  .\Ellissi\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
            #        }
            #        if (Test-Path ".\graph_Parametri_Ellissi_Funz_Ottiche_T.png")
	        #        {
            #           Move-Item -Force graph_Parametri_Ellissi_Funz_Ottiche_T.png  .\Ellissi_T\Def_${grad_d}__Foc_${grad_f}__Lelem_${lung_elem}__Ldri_${lung_drift_m}\Energia_${energy}_MeV.png
            #        }                    
                    
            #        $g=$g+1
            #    }
            #    $energy=${ENERGY}
            #}