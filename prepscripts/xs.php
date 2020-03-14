<?php

$dir = "/home/ben/activity/xs/tendl";
$out_dir = "/home/ben/activity/xs/datatemp";
$root = scandir($dir); 




//Read files

foreach($root as $value){ 
	if(stristr($value,".tendl")){
		
		//Process tendl file
		process_file($dir, $value, $out_dir);

	}
}



function process_file($file_dir, $file_name, $out_dir){

	$file_path = $file_dir."/".$file_name;
	echo "Processing: ".$file_path.chr(10);
	
	$output_file = "";
	
	$cross_section_points = 0;
	$product = 0;
	
	//projectile
	if(stristr($file_name,"-p")){
		$projectile_z = 1;
		$projectile_a = 1;
	}

	$fh = fopen($file_path, 'r');
    while(!feof($fh)){
		$file_row = fgets($fh);
        $file_row = str_replace(chr(10),"",$file_row);
        $file_row = str_replace(chr(10),"",$file_row);
		//echo $file_row.chr(10);
		
		$mf = 1 * substr($file_row,70,2);
		$mt = 1 * substr($file_row,72,3);
		$row = 1 * substr($file_row,75,5);
		
		//read in target
		if($mf==1&&$mt==451&&$row==1){
			$data_array = read_data_line($file_row, true);
			$target = $data_array[1];
			$target_z = floor($data_array[1]/1000);
			$target_a = $target - 1000 * $target_z;
			if(stristr($file_name,"m-")){
				$target_m = 1;
			}else{
				$target_m = 0;
			}
			
			$output_file_name = $projectile_z."_".$projectile_a."-".$target_z."_".$target_a."_".$target_m.".dat";
		}
		
		//Read in cross section data
		if($mf==3&&$mt==5&&$row==3){
			$data_array = read_data_line($file_row);
			$cross_section_points = $data_array[1];
		}
		
		//Read in cross section data
		if($cross_section_points>0){
			$rows_to_read = ceil($cross_section_points/3);
			$data_counter = 0;
			for($i=1; $i<=$rows_to_read; $i++){
				$file_row = fgets($fh);
				$data_array = read_data_line($file_row, true);
				for($j=1; $j<=3; $j++){
					$data_counter = $data_counter + 1;
					if($data_counter<=$cross_section_points){
						$cross_section[$data_counter]->x = $data_array[2*$j-1];
						$cross_section[$data_counter]->y = $data_array[2*$j];
					}
				}				
			}
			//end
			$cross_section_points = 0;
		}
		
		//product yield
		if($mf==6&&$mt==5){
			$data_array = read_data_line($file_row, true);
			if($data_array[3]==0&&$data_array[4]==1&&$data_array[5]==1){
				$product = $data_array[1];
			}
		}
		
		//product yield
		if($product>0){
			$product_z = floor($product/1000);
			$product_a = $product - 1000 * $product_z;
			$product_m = 0;
			//echo $z." ".$a." ".$m." ".chr(10);
			
			$output_file .= "#Header ".chr(10);
			$output_file .= "#Target ".$target_z." ".$target_a." ".$target_m.chr(10);
			$output_file .= "#Product ".$product_z." ".$product_a." ".$product_m.chr(10);
			
			$file_row = fgets($fh);
			$data_array = read_data_line($file_row, true);	
			$data_points = $data_array[1];
			$rows_to_read = ceil($data_points/3);
			$output_file .= "#Datapoints ".$data_points.chr(10);	//store data point count
			$data_counter = 0;
			for($i=1; $i<=$rows_to_read; $i++){
				$file_row = fgets($fh);
				$data_array = read_data_line($file_row, true);
				for($j=1; $j<=3; $j++){
					$data_counter = $data_counter + 1;
					if($data_counter<=$data_points){
						//store data temporarily
						$data_x = $data_array[2*$j-1];
						$data_y = $data_array[2*$j];
						//multiply by total cross section
						$set = false;
						for($k=1; $k<=count($cross_section); $k++){
							if($data_x==$cross_section[$k]->x){
								$data_y = $data_y * $cross_section[$k]->y;
								$set = true;
							}
						}
						if($set==false){
							if($data_x<$cross_section[1]->x||$data_x>$cross_section[count($cross_section)]->x){
								$data_y = 0;
							}else{
								for($k=2; $k<=count($cross_section); $k++){
									if($data_x>$cross_section[$k-1]->x&&$data_x<$cross_section[$k]->x){
										$data_y = $data_y * $cross_section[$k]->y;
										$set = true;
									}
								}
							}
						}
						//store data
						$output_file .= $data_x." ".$data_y.chr(10);
					}
				}
			}
			
			//end
			$product = 0;
		}
	}
	fclose($fh);
	
	//echo $output_file;
	//print_r($cross_section);

	$fh = fopen($out_dir."/".$output_file_name, 'w');
	fwrite($fh, $output_file);
	fclose($fh);
		
}

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