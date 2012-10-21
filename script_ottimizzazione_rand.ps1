$GRADIENTI_INIZIALI=50.0
$GRADIENTI_FINALI=150.0
$NUM_CICLI=100000
$L_TOT=0.45
$L_IN=0.2
$L_ELEM=0.06
$ENERGY=30
$DRIFT_F=0.8

if (!(Test-Path "Plot")){
new-Item -type directory ("Plot")}

if (!(Test-Path "Input")){
new-Item -type directory ("Input")}

$a=0

while ($a -le $NUM_CICLI)
{
                
                $L_1= Get-Random -min 0.001 -max 0.1
                $L_1=[Math]::Round($L_1, 3)
                $L_2= Get-Random -min 0.001 -max 0.05
                $L_2=[Math]::Round($L_2, 3)
                $grad_1= Get-Random -min 60 -max 140
                $grad_1=[Math]::Round($grad_1, 1)
                $grad_2= Get-Random -min 60 -max 140
                $grad_2=[Math]::Round($grad_2, 1)
                $grad_3= Get-Random -min 60 -max 140
                $grad_3=[Math]::Round($grad_3, 1)            
                
                
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
                $a=$a+1
                echo $a
}
