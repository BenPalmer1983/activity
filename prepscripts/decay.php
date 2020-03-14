<?php

$dir = "/home/ben/activity/data/decay";
$out_dir = "/home/ben/activity/data";
$root = scandir($dir); 

//preparation - half lives
$i = 0;
$half_life_unit[$i]->unit = "KEV";
$half_life_unit[$i]->factor = 1.77E-25;
$i++;
$half_life_unit[$i]->unit = "MEV";
$half_life_unit[$i]->factor = 1.77E-22;
$i++;
$half_life_unit[$i]->unit = "NS";
$half_life_unit[$i]->factor = 1E-9;
$i++;
$half_life_unit[$i]->unit = "MS";
$half_life_unit[$i]->factor = 0.001;
$i++;
$half_life_unit[$i]->unit = "S";
$half_life_unit[$i]->factor = 1;
$i++;
$half_life_unit[$i]->unit = "M";
$half_life_unit[$i]->factor = 60;
$i++;
$half_life_unit[$i]->unit = "H";
$half_life_unit[$i]->factor = 3600;
$i++;
$half_life_unit[$i]->unit = "D";
$half_life_unit[$i]->factor = 1*24*3600;
$i++;
$half_life_unit[$i]->unit = "Y";
$half_life_unit[$i]->factor = 365*1*24*3600;

for($i=0; $i<count($half_life_unit); $i++){
	$half_life_unit[$i]->unit = strtoupper(" ".$half_life_unit[$i]->unit." ");
}


//preparation - decay type
$i = 0;
$decay_mode_array[$i]->unit = "A";	//Alpha
$decay_mode_array[$i]->dZ = -2;
$decay_mode_array[$i]->dA = -4;
$i++;
$decay_mode_array[$i]->unit = "N";	//Neutron
$decay_mode_array[$i]->dZ = 0;
$decay_mode_array[$i]->dA = -1;
$i++;
$decay_mode_array[$i]->unit = "B-";	//Electron
$decay_mode_array[$i]->dZ = 1;
$decay_mode_array[$i]->dA = 0;
$i++;
$decay_mode_array[$i]->unit = "B+";	//Electron
$decay_mode_array[$i]->dZ = -1;
$decay_mode_array[$i]->dA = 0;
$i++;
$decay_mode_array[$i]->unit = "BN";	//Electron + Neutron
$decay_mode_array[$i]->dZ = 1;
$decay_mode_array[$i]->dA = -1;
$i++;
$decay_mode_array[$i]->unit = "EC+B+";	//Electron Capture + Positron Emission
$decay_mode_array[$i]->dZ = -1;
$decay_mode_array[$i]->dA = 0;
$i++;
$decay_mode_array[$i]->unit = "EP";	//Positron Emission + Proton
$decay_mode_array[$i]->dZ = -2;
$decay_mode_array[$i]->dA = -1;
$i++;
$decay_mode_array[$i]->unit = "EC";	//Electron Capture
$decay_mode_array[$i]->dZ = -1;
$decay_mode_array[$i]->dA = 0;
$i++;
$decay_mode_array[$i]->unit = "EP AP";	//Positron Emission + Proton
$decay_mode_array[$i]->dZ = -2;
$decay_mode_array[$i]->dA = -1;
$i++;
$decay_mode_array[$i]->unit = "P";	//Proton
$decay_mode_array[$i]->dZ = -1;
$decay_mode_array[$i]->dA = -1;
$i++;
$decay_mode_array[$i]->unit = "2P";	//Proton
$decay_mode_array[$i]->dZ = -2;
$decay_mode_array[$i]->dA = -2;



$stable = "";
$unstable = "";
$gamma_lines = "";




//Read files

foreach($root as $value){ 

	//if(stristr($value,"decay_0540_27-Co-55.dat")){
	if(stristr($value,".dat")){
		//echo "Processing: ".$value.chr(10);
		
		$file_name = str_replace(".dat","",$value);
		$file_name_array = explode("_",$file_name);
		$atom_name = explode('-',$file_name_array[2]);
			
		$z = $atom_name[0];
		$element = $atom_name[1];
		$a = $atom_name[2];
		$meta = 0;
		if(stristr($a,"M")){
			$meta = 1;
			$a = str_replace("M","",$a);
		}
		if(stristr($a,"N")){
			$meta = 2;
			$a = str_replace("N","",$a);
		}
			
		//echo $z.",".$element.",".$a.",".$meta.chr(10);
		
		$fh = fopen($dir."/".$value, 'r');
		$file_data = fread($fh, filesize($dir."/".$value));
		
		$file_data = explode(chr(10),$file_data);
		for($i=0; $i<count($file_data); $i++){
			$file_row = $file_data[$i];
			$file_row_data = substr($file_row,0,66);
			//Half Life
			if(stristr($file_row_data,"Parent half-life:")){
			
				unset($temp_array);
				$temp_array = explode("Parent half-life:",$file_row_data);
				$half_life = strtoupper(trim($temp_array[1]));
				
				$half_life_value = 0;
				if($half_life=="STABLE"||$half_life=="[STABLE]"){
					$half_life_value = -1;
					$stable .= $a." ".$z." ".$meta.chr(10);
				}else{
					//echo "[".$half_life."] ";
					for($j=0; $j<count($half_life_unit); $j++){
						if(stristr($half_life,$half_life_unit[$j]->unit)){
							$temp_half_life_array = explode($half_life_unit[$j]->unit,$half_life);
							//echo "(_".$temp_half_life_array[0].")".chr(10);
							$half_life_value = $temp_half_life_array[0] * $half_life_unit[$j]->factor;
							$j = count($half_life_unit);
						}
					}
				}
				
				//echo $half_life_value.chr(10);
				
				
			}
			//Decay Modes
			if(stristr($file_row_data,"Decay Mode:")){
				$temp_array = explode("Decay Mode:",$file_row_data);
				$decay_type = trim($temp_array[1]);
				$decay_array = explode(",",$decay_type);
				for($j=0;$j<count($decay_array);$j++){
					$decay_mode = str_replace("%","",trim($decay_array[$j]));
					if(stristr($decay_mode,"=")){
						$decay_mode_temp = explode("=",$decay_mode);
						for($k=0; $k<count($decay_mode_array); $k++){
							if($decay_mode_temp[0]==$decay_mode_array[$k]->unit){
							
								$decay_probability = explode(" ",$decay_mode_temp[1]);
								$decay_probability = $decay_probability[0] * 1;
							
								//echo $decay_mode.",".$half_life_value.",".$decay_mode_temp[0].chr(10);
								//echo $a." ".$z." ".$meta." ".($a+$decay_mode_array[$k]->dA)." ".($z+$decay_mode_array[$k]->dZ)." ".($decay_probability/100)." ".$half_life_value." ".$decay_mode_temp[0].chr(10);
							
								if($half_life_value>=0){
									$unstable .= strtoupper($element)." ".$a." ".$z." ".$meta." ".($a+$decay_mode_array[$k]->dA)." ".($z+$decay_mode_array[$k]->dZ)." ".($decay_probability/100)." ".$half_life_value." ".$decay_mode_temp[0].chr(10);
								}
							}
						}
					}
				}
				/*
				$half_life_unit[$i]->unit = "EP AP";	//Positron Emission + Proton
				$half_life_unit[$i]->dZ = -2;
				$half_life_unit[$i]->dA = -1;
				*/
				//decay_mode_array
				//echo $file_row.chr(10);
			}
			
			//echo $file_data[$i].chr(10);
			
			
			//gamma lines
			$mf = 1 * substr($file_row,70,2);
			$mt = 1 * substr($file_row,72,3);
			$row = 1 * substr($file_row,75,5);
			
			if($mf==8&&$mt==457){				
				$row_data = read_data_line($file_row, true);
				if($row_data[1]==0&&$row_data[2]==0&&$row_data[3]==0&&$row_data[4]==0&&$row_data[5]>=6&&$row_data[6]>=1&&$row_data[5]<=30&&$row_data[6]<=10000){
				if((($row_data[5]/6)-round($row_data[5]/6))==0){
					
					$data_field_count = $row_data[6];	
					
					$gamma_lines .= "#Header ".chr(10);
					$gamma_lines .= "#Parent ".$z." ".$a." ".$meta.chr(10);
					$gamma_lines .= "#Datapoints ".$data_field_count.chr(10);
										
					//skip data rows
					for($j=1; $j<=($row_data[5]/6); $j++){
						$i++;
					}
					
					for($j=1; $j<=$data_field_count; $j++){
						$i++;
						//data header
						$file_row = $file_data[$i];
						$row_data = read_data_line($file_row, true);
						$rows_of_data = $row_data[5];
						$gamma_energy = $row_data[1];
												
						for($k=1; $k<=($rows_of_data/6); $k++){							
							$i++;
							$file_row = $file_data[$i];
							$row_data = read_data_line($file_row, true);
							if($k==1){
								$gamma_intensity = $row_data[3];
							}
						}
						
						$gamma_lines .= $gamma_energy." ".$gamma_intensity.chr(10);
					
					}
					
					
					
					
					
					
					
				}	
				}
			}
			
		}
		
	}
}



$output_file_name = "decaymodes.txt";

//$fh = fopen($out_dir."/".$output_file_name, 'w');
//fwrite($fh, $unstable);
//fclose($fh);

$output_file_name = "gamma_energies.txt";

$fh = fopen($out_dir."/".$output_file_name, 'w');
fwrite($fh, $gamma_lines);
fclose($fh);









function read_data_line($input, $format=false){	

	if($format==true){
		$input = str_replace(chr(13).chr(10),"",$input);
		$input = str_replace(chr(13),"",$input);
		$input = str_replace(chr(10),"",$input);
	}		
	$data_array[1] = substr($input,0,11);
	$data_array[2] = substr($input,11,11);
	$data_array[3] = substr($input,22,11);
	$data_array[4] = substr($input,33,11);
	$data_array[5] = substr($input,44,11);
	$data_array[6] = substr($input,55,11);
	if($format==true){
	
		for($i=1; $i<=6; $i++){
			$factor = 1;
			if(substr($data_array[$i],0,1)=="-"){
				$factor = -1;
				$data_array[$i] = substr($data_array[$i],1,strlen($data_array[$i])-1);
			}
			$data_array[$i] = str_replace("-","E-",$data_array[$i]);
			$data_array[$i] = str_replace("+","E+",$data_array[$i]);
			$data_array[$i] = $factor * $data_array[$i];
		}
	}	
	return $data_array;
}





?>