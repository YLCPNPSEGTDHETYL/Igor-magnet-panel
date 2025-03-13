Menu "Load Waves"
	"Load Hall Data File...", LoadHalldata()
End

function Load_rampfield_data(dataname)
	string dataname
	
	newpath/O/Q pathname
	
	if (V_flag != 0)
	print "Loading canceled."
		return 0
	endif
		
	pathinfo pathname
	string pathnamebase=S_path
	
	variable index=0
	
	string filename
    string kValue,kvalue2
    
    //variable startIndex
    //variable endIndex,endIndex2
    
    variable kpos,InitialIndex
    
    variable fileCount = 0
    variable LoadCount = 0
    
    Do
    	filename = indexedfile(pathname, index, "????")
		if (strlen(filename) == 0)
			break
		endif
		if (strsearch(filename, "K", 0) >= 0 && strsearch(filename, dataname, 0) >= 0)
			Print "Loading file..."
			LoadCount += 1
			break
		else
		endif
		index += 1
    while(1)
    
    if (LoadCount == 0)
    	Print "Loading file..."
    	Print "Could not find any files with "+ dataname +" in the name"
    	return -1
    endif    
    
    newdatafolder/O/S $dataname
    
    
    make/O/D/N=0 kValues
	
	do
        filename = indexedfile(pathname, index, "????")
        if (strlen(filename) == 0)
            break
        endif

        if (strsearch(filename, "K", 0) >= 0 && strsearch(filename, dataname, 0) >= 0)
            //startIndex = strsearch(filename, "_", 0) + 1
            //endIndex = strsearch(filename, "K", startIndex) + 1
            //kValue = filename[startIndex, endIndex - 1]   
            
            //endIndex2 = strsearch(filename, "K", startIndex)
            //kValue2 = filename[startIndex, endIndex2 - 1]
            
            kpos = strsearch(filename, "K_ramp", 0)
            InitialIndex = strsearch(filename, "_", kpos, 1)
            
            kValue = filename[InitialIndex + 1, kpos]
            kValue2 = filename[InitialIndex + 1, kpos - 1]
            
            newdatafolder/O/S $kValue
            
            loadwave/A/Q/G/D/P=pathname filename
            
            
            wave wave0, wave1, wave2, wave3, wave4, wave5
            
            Duplicate/O wave0, $(kValue + "_time_stamp_rawdata"); KillWaves wave0
            Duplicate/O wave1, $(kValue + "_Temp_chA_rawdata"); KillWaves wave1
            Duplicate/O wave2, $(kValue + "_Temp_chB_rawdata"); KillWaves wave2
            Duplicate/O wave3, $(kValue + "_field_rawdata"); KillWaves wave3
            Duplicate/O wave4, $(kValue + "_voltage_rawdata"); KillWaves wave4
            Duplicate/O wave5, $(kValue + "_resistance_rawdata"); KillWaves wave5
            
            Note $(kValue + "_time_stamp_rawdata"), "time (s)"
            Note $(kValue + "_Temp_chA_rawdata"), "Temperature (K)"
            Note $(kValue + "_Temp_chB_rawdata"), "Temperature (K)"
            Note $(kValue + "_field_rawdata"), "Magnetic field (T)"
            Note $(kValue + "_voltage_rawdata"), "Voltage (V)"
            Note $(kValue + "_resistance_rawdata"), "Resistance (Ohm)"
            
            
			wave fieldwave = $(kValue + "_field_rawdata")
			wave resistancewave = $(kValue + "_resistance_rawdata")
			
			
			ControlInfo CutValue
			variable CutValue = V_Value
			
			if (CutValue == 1)
            	CutZerodata(kValue, fieldwave)
            endif
            
            ControlInfo CalcRyxValue
            variable CalcRyxValue = V_Value
            
            if (CalcRyxValue == 1)
            	calc_Ryx(kValue,fieldwave,resistancewave)
            endif
            
            setdatafolder ::
            
                        
            InsertPoints (fileCount),1, kValues
            kValues[fileCount] = str2num(kValue2)
                        
            fileCount += 1
        endif        
        index += 1
    while(1)
    
    Sort kValues,kValues
    
    setdatafolder ::
	if (CutValue == 1)
		Print "Extra zero point removed."
	endif
	if (CalcRyxValue == 1)
		Print "Rxx & Ryx Calculations are completed."
	endif
    
end

function calc_Ryx(kValue,fieldwave,resistancewave)
	string kValue
	wave fieldwave
	wave resistancewave		
	
	CountSplit(kValue,fieldwave,resistancewave)
		
	//up
	
	Duplicate/O $(kValue + "_field_pu"), up_field, Ryx_u, Rxx_u
	Duplicate/O $(kValue + "_resistance_pu"), resistance_pu_dummy
	Duplicate/O $(kValue + "_field_mu"), field_mu_dummy
	field_mu_dummy = -field_mu_dummy
	
	Ryx_u = ( resistance_pu_dummy - interp(up_field, field_mu_dummy, $(kValue + "_resistance_mu")) ) / 2
	Rxx_u = ( resistance_pu_dummy + interp(up_field, field_mu_dummy, $(kValue + "_resistance_mu")) ) / 2
	
	killwaves field_mu_dummy, resistance_pu_dummy
	
	
	//down
	
	Duplicate/O $(kValue + "_field_pd"), down_field, Ryx_d, Rxx_d
	Duplicate/O $(kValue + "_resistance_pd"), resistance_pd_dummy
	Duplicate/O $(kValue + "_field_md"), field_md_dummy
	field_md_dummy = -field_md_dummy
	
	Ryx_d = ( resistance_pd_dummy - interp(down_field, field_md_dummy, $(kValue + "_resistance_md")) ) / 2
	Rxx_d = ( resistance_pd_dummy + interp(down_field, field_md_dummy, $(kValue + "_resistance_md")) ) / 2
	
	killwaves field_md_dummy, resistance_pd_dummy
	
	make/D/O/N=(numpnts(up_field)+numpnts(down_field)) fieldall,Ryx, Rxx
	
	variable ii = 0	
	do
		fieldall[ii] = -down_field[ii]
		Ryx[ii] = -Ryx_d[ii]
		Rxx[ii] = Rxx_d[ii]
		ii += 1
	while(ii < numpnts(down_field))
	
	ii = 0
	do
		fieldall[ii + numpnts(down_field)] = up_field[ii]
		Ryx[ii + numpnts(down_field)] = Ryx_u[ii]
		Rxx[ii + numpnts(down_field)] = Rxx_u[ii]
		ii += 1
	while(ii < numpnts(up_field))
	
	Note/K Ryx; Note/K Rxx; Note/K fieldall;
	Note Ryx, "Resistance (Ohm)"; Note Rxx, "Resistance (Ohm)"; Note fieldall, "Magnetic field (T)"
	
	Duplicate/O Ryx, $(kValue + "_Ryx"); KillWaves Ryx
    Duplicate/O Rxx, $(kValue + "_Rxx"); KillWaves Rxx
    Duplicate/O fieldall, $(kValue + "_field"); KillWaves fieldall
	
	killwaves up_field,down_field,Rxx_d,Rxx_u,Ryx_u,Ryx_d
	killwaves $(kValue + "_field_pu"),$(kValue + "_field_pd"),$(kValue + "_field_mu"),$(kValue + "_field_md")
	killwaves $(kValue + "_resistance_pu"),$(kValue + "_resistance_pd"),$(kValue + "_resistance_mu"),$(kValue + "_resistance_md")
	
end

function calc_Ryx_kValues()
	
	wave kValues
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
		
	do			
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue		
		
		wave fieldwave = $(kValue + "_field_rawdata")
		wave resistancewave = $(kValue + "_resistance_rawdata")
		
		CountSplit(kValue,fieldwave,resistancewave)
			
		//up
			
		Duplicate/O $(kValue + "_field_pu"), up_field, Ryx_u, Rxx_u
		Duplicate/O $(kValue + "_resistance_pu"), resistance_pu_dummy
		Duplicate/O $(kValue + "_field_mu"), field_mu_dummy
		field_mu_dummy = -field_mu_dummy
		
		Ryx_u = ( resistance_pu_dummy - interp(up_field, field_mu_dummy, $(kValue + "_resistance_mu")) ) / 2
		Rxx_u = ( resistance_pu_dummy + interp(up_field, field_mu_dummy, $(kValue + "_resistance_mu")) ) / 2
		
		killwaves field_mu_dummy, resistance_pu_dummy
		
		
		//down
		
		Duplicate/O $(kValue + "_field_pd"), down_field, Ryx_d, Rxx_d
		Duplicate/O $(kValue + "_resistance_pd"), resistance_pd_dummy
		Duplicate/O $(kValue + "_field_md"), field_md_dummy
		field_md_dummy = -field_md_dummy
		
		Ryx_d = ( resistance_pd_dummy - interp(down_field, field_md_dummy, $(kValue + "_resistance_md")) ) / 2
		Rxx_d = ( resistance_pd_dummy + interp(down_field, field_md_dummy, $(kValue + "_resistance_md")) ) / 2
		
		killwaves field_md_dummy, resistance_pd_dummy
		
		make/D/O/N=(numpnts(up_field)+numpnts(down_field)) fieldall,Ryx, Rxx
		
		variable ii = 0	
		do
			fieldall[ii] = -down_field[ii]
			Ryx[ii] = -Ryx_d[ii]
			Rxx[ii] = Rxx_d[ii]
			ii += 1
		while(ii < numpnts(down_field))
		
		ii = 0
		do
			fieldall[ii + numpnts(down_field)] = up_field[ii]
			Ryx[ii + numpnts(down_field)] = Ryx_u[ii]
			Rxx[ii + numpnts(down_field)] = Rxx_u[ii]
			ii += 1
		while(ii < numpnts(up_field))
		
		Note/K Ryx; Note/K Rxx; Note/K fieldall;
		Note Ryx, "Resistance (Ohm)"; Note Rxx, "Resistance (Ohm)"; Note fieldall, "Magnetic field (T)"
	
		Duplicate/O Ryx, $(kValue + "_Ryx"); KillWaves Ryx
    	Duplicate/O Rxx, $(kValue + "_Rxx"); KillWaves Rxx
    	Duplicate/O fieldall, $(kValue + "_field"); KillWaves fieldall
		
		killwaves up_field,down_field,Rxx_d,Rxx_u,Ryx_u,Ryx_d
		killwaves $(kValue + "_field_pu"),$(kValue + "_field_pd"),$(kValue + "_field_mu"),$(kValue + "_field_md")
		killwaves $(kValue + "_resistance_pu"),$(kValue + "_resistance_pd"),$(kValue + "_resistance_mu"),$(kValue + "_resistance_md")	
		
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	Print "Rxx & Ryx Calculations are completed."
end


function calc_rho_yx_kValues()
	
	wave kValues
	
	String currentDF = GetDataFolder(1)
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
	
	do			
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue
		
		ControlInfo ThkVar_Alltab
		variable ThkVar = V_value
		
		
		if (ThkVar == 0)
			SetDataFolder currentDF
			Print "Error: Received an invalid input. A non-zero finite thickness value must be entered."
			Abort "Error: \nReceived an invalid input. A non-zero finite thickness value must be entered."
		return 0
		endif
		
		if (!WaveExists($(kValue + "_Ryx")))
			Print "Error: Expected Ryx wave, but it does not exist. Please calculate Ryx first."
			SetDataFolder currentDF
		return 0
		endif
		
		wave Ryx = $(kValue + "_Ryx")
		Duplicate/O Ryx, rho_yx
		
		rho_yx = Ryx * ThkVar * 10^(-1) // mm --> cm
		
		Note/K rho_yx
		Note rho_yx, "Resistivity (Ohm cm)"
		Note rho_yx, "Thickness: t = " + num2str(ThkVar) + " mm"
		
		Duplicate/O rho_yx, $(kValue + "_rho_yx"); KillWaves rho_yx
		
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	
	Print "rho_yx Calculations are completed."
end


function calc_RH()
	
	wave kValues
	
	String currentDF = GetDataFolder(1)
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
		
	do			
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue		
		
		wave fieldwave = $(kValue + "_field")
		wave rhoyxwave = $(kValue + "_rho_yx")
		
		
		Duplicate/O fieldwave, RH
		
		RH = rhoyxwave/fieldwave
		
		Note/K RH
		Note RH, "Hall coefficient (cm^3/C) = (Ohm cm/T)"
		
		
    	Duplicate/O RH, $(kValue + "_RH"); KillWaves RH
    	
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	Print "RH Calculations are completed."
end


function calc_carrier_density()
	
	wave kValues
	
	String currentDF = GetDataFolder(1)
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
		
	variable echaerge = -1.60217663 * 10^(-19) /// C
	do			
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue		
		
		wave RHwave = $(kValue + "_RH")
		
		
		Duplicate/O RHwave, carrier_density
		
		
		carrier_density = 1/ (RHwave * echaerge)
		
		Note/K carrier_density
		Note carrier_density, "Carrier density (1/cm^3)"
		
		
    	Duplicate/O carrier_density, $(kValue + "_carrier_density"); KillWaves carrier_density
    	
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	Print "carrier density Calculations are completed."
end


function calc_rho_xx_kValues()
	
	wave kValues
	
	String currentDF = GetDataFolder(1)
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
		
	do			
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue		
		
		wave fieldwave = $(kValue + "_field_rawdata")
		wave resistancewave = $(kValue + "_resistance_rawdata")
		
		CountSplit(kValue,fieldwave,resistancewave)
			
		Sort $(kValue + "_field_pd"),$(kValue + "_field_pd"),$(kValue + "_resistance_pd")
		Sort $(kValue + "_field_mu"),$(kValue + "_field_mu"),$(kValue + "_resistance_mu")
			
		Duplicate/O $(kValue + "_field_pu"), plus_field
		Duplicate/O $(kValue + "_field_mu"), minus_field		
		Duplicate/O $(kValue + "_resistance_pu"), plus_resistance, resistance_pu
		Duplicate/O $(kValue + "_resistance_mu"), minus_resistance, resistance_mu
		
		
		
		plus_resistance = ( resistance_pu + interp(plus_field, $(kValue + "_field_pd"), $(kValue + "_resistance_pd")) )/2
		minus_resistance = ( resistance_mu + interp(minus_field, $(kValue + "_field_md"), $(kValue + "_resistance_md")) )/2
		
		
		make/D/O/N=(numpnts(plus_field)+numpnts(minus_field)) fieldall,resistanceall
		
		variable ii = 0	
		do
			fieldall[ii] = minus_field[ii]
			resistanceall[ii] = minus_resistance[ii]
			ii += 1
		while(ii < numpnts(minus_field))
		
		ii = 0
		do
			fieldall[ii + numpnts(minus_field)] = plus_field[ii]
			resistanceall[ii + numpnts(minus_field)] = plus_resistance[ii]
			ii += 1
		while(ii < numpnts(plus_field))
		
		ControlInfo ThkVar_Alltab
		variable ThkVar = V_value
		
		ControlInfo LenVar_DAtab1
		variable LenVar = V_value
		
		ControlInfo WidVar_DAtab1
		variable WidVar = V_value
		
		if (ThkVar == 0 || LenVar == 0 || WidVar == 0)
			SetDataFolder currentDF
			Print "Error: Received an invalid input. A non-zero finite sample size must be entered."
			Abort "Error: \nReceived an invalid input. A non-zero finite sample size must be entered."
		return 0
		endif
		
		Duplicate/O fieldall, rho_xx
		
		rho_xx = (resistanceall * ThkVar * WidVar) /LenVar * 10^(-1) // mm --> cm
		
		Note/K rho_xx
		Note rho_xx, "Resistivity (Ohm cm)"
		Note rho_xx, "Thickness: t = " + num2str(ThkVar) + " mm"
		Note rho_xx, "Length: L = " + num2str(LenVar) + " mm"
		Note rho_xx, "Width: w = " + num2str(WidVar) + " mm"
		
		
		
		Note/K resistanceall; Note/K fieldall; Note/K rho_xx;
		Note resistanceall, "Resistance (Ohm)"; Note fieldall, "Magnetic field (T)"; Note rho_xx, "Resistivity (Ohm cm)"
		
    	Duplicate/O resistanceall, $(kValue + "_resistance"); KillWaves resistanceall
    	Duplicate/O fieldall, $(kValue + "_field"); KillWaves fieldall
    	Duplicate/O rho_xx, $(kValue + "_rho_xx"); KillWaves rho_xx
		
		killwaves plus_field,minus_field,plus_resistance,minus_resistance
		killwaves resistance_mu,resistance_pu
		killwaves $(kValue + "_field_pu"),$(kValue + "_field_pd"),$(kValue + "_field_mu"),$(kValue + "_field_md")
		killwaves $(kValue + "_resistance_pu"),$(kValue + "_resistance_pd"),$(kValue + "_resistance_mu"),$(kValue + "_resistance_md")	
		
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	Print "rho_xx Calculations are completed."
end


function calc_MR_kValues()
	
	ControlInfo Remove_yx_HDtab1
	Variable Remove_yx = V_value
		
	if (Remove_yx == 1)
		calc_Ryx_kValues()
	endif
	
	wave kValues
	
	String currentDF = GetDataFolder(1)
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
		
	do			
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue			
		
		wave fieldwave = $(kValue + "_field_rawdata")
		wave resistancewave = $(kValue + "_resistance_rawdata")
		
		CountSplit(kValue,fieldwave,resistancewave)
			
		Sort $(kValue + "_field_pd"),$(kValue + "_field_pd"),$(kValue + "_resistance_pd")
		Sort $(kValue + "_field_mu"),$(kValue + "_field_mu"),$(kValue + "_resistance_mu")
			
		Duplicate/O $(kValue + "_field_pu"), plus_field
		Duplicate/O $(kValue + "_field_mu"), minus_field		
		Duplicate/O $(kValue + "_resistance_pu"), plus_resistance, resistance_pu
		Duplicate/O $(kValue + "_resistance_mu"), minus_resistance, resistance_mu
		
		
		
		plus_resistance = ( resistance_pu + interp(plus_field, $(kValue + "_field_pd"), $(kValue + "_resistance_pd")) )/2
		minus_resistance = ( resistance_mu + interp(minus_field, $(kValue + "_field_md"), $(kValue + "_resistance_md")) )/2
		
		
		
		make/D/O/N=(numpnts(plus_field)+numpnts(minus_field)) fieldall,resistanceall
		
		variable ii = 0	
		do
			fieldall[ii] = minus_field[ii]
			resistanceall[ii] = minus_resistance[ii]
			ii += 1
		while(ii < numpnts(minus_field))
		
		ii = 0
		do
			fieldall[ii + numpnts(minus_field)] = plus_field[ii]
			resistanceall[ii + numpnts(minus_field)] = plus_resistance[ii]
			ii += 1
		while(ii < numpnts(plus_field))
		
		
		Duplicate/O fieldall, MR
		
		
		Note/K fieldall; Note/K MR;
		Note fieldall, "Magnetic field (T)"
		
		MR = ( resistanceall - interp(0,fieldall,resistanceall) ) /(interp(0,fieldall,resistanceall))
		
		
    	Duplicate/O resistanceall, $(kValue + "_resistance"); KillWaves resistanceall
    	Duplicate/O fieldall, $(kValue + "_field"); KillWaves fieldall
		
		killwaves plus_field,minus_field,plus_resistance,minus_resistance
		killwaves resistance_mu,resistance_pu
		killwaves $(kValue + "_field_pu"),$(kValue + "_field_pd"),$(kValue + "_field_mu"),$(kValue + "_field_md")
		killwaves $(kValue + "_resistance_pu"),$(kValue + "_resistance_pd"),$(kValue + "_resistance_mu"),$(kValue + "_resistance_md")	
		
		
		if (Remove_yx == 1)
			Killwaves $(kValue + "_resistance")
			Duplicate/O $(kValue + "_Rxx"),Rxx
			
			MR = ( Rxx - interp(0,$(kValue + "_field"),Rxx) ) /(interp(0,$(kValue + "_field"),Rxx))
				
			Killwaves Rxx
		endif
		
    	Duplicate/O MR, $(kValue + "_MR"); KillWaves MR
		
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	Print "MR Calculations are completed."
end

function calc_field_square()
	
	wave kValues
	
	String currentDF = GetDataFolder(1)
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
		
	do			
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue		
		
		wave fieldwave = $(kValue + "_field")
		
		
		Duplicate/O fieldwave, field_square
		
		field_square = fieldwave*fieldwave
		
		Note/K field_square
		Note field_square, "sqrt(Magnetic field) (T^2)"
		
		
    	Duplicate/O field_square, $(kValue + "_field_square"); KillWaves field_square
    	
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	Print "field_square Calculations are completed."
end


function CountSplit(kValue,fieldwave,resistancewave)
	string kValue	
	wave fieldwave
	wave resistancewave
	
	if (fieldwave[0] == 0)
        CutZerodata(kValue, fieldwave)
    endif
	
	variable numPoints = numpnts(fieldwave)
	variable ii = 0
	
	variable segmentState = 1
	variable firstCount = 0, secondCount = 0, thirdCount = 0, fourthCount = 0
	
	do   
        if (segmentState == 1)
        
            if (fieldwave[ii] < fieldwave[ii+1] && 0 <= fieldwave[ii+1])
                //firstCount += 1
            else
                segmentState = 2
                firstCount = ii-1
                //secondCount = firstCount
            endif
            
            
        elseif (segmentState == 2)
        
            if (fieldwave[ii] > fieldwave[ii+1] && 0 <= fieldwave[ii])
                //secondCount += 1
            else
                segmentState = 3
                //thirdCount = secondCount
                //thirdCount += 1
                secondCount = ii-1
            endif
            
            
        elseif (segmentState == 3)
        
            if (fieldwave[ii] > fieldwave[ii+1] && 0 > fieldwave[ii])
                //thirdCount += 1
            else
                segmentState = 4
                //thirdCount += 1
                //fourthCount = thirdCount
                thirdCount = ii-1
            endif
            
            
        elseif (segmentState == 4)
        
            if (fieldwave[ii-1] < fieldwave[ii] && 0 > fieldwave[ii])
                //fourthCount += 1
            endif
        endif


        ii += 1
        fourthCount = ii-1
    while(ii < numPoints)
	
	// plus up
	Duplicate/O/R=[0,(firstCount)] fieldwave, $(kValue + "_field_pu")
	Duplicate/O/R=[0,(firstCount)] resistancewave, $(kValue + "_resistance_pu")
	// plus down
	Duplicate/O/R=[(firstCount+1),(secondCount)] fieldwave, $(kValue + "_field_pd")
	Duplicate/O/R=[(firstCount+1),(secondCount)] resistancewave, $(kValue + "_resistance_pd")
	// minus up
	Duplicate/O/R=[(secondCount+1),(thirdCount)] fieldwave, $(kValue + "_field_mu")
	Duplicate/O/R=[(secondCount+1),(thirdCount)] resistancewave, $(kValue + "_resistance_mu")
	// minus down
	Duplicate/O/R=[(thirdCount+1),(fourthCount)] fieldwave, $(kValue + "_field_md")
	Duplicate/O/R=[(thirdCount+1),(fourthCount)] resistancewave, $(kValue + "_resistance_md")
	
	
end


function CutZerodata(kValue,fieldwave)
	string kValue
	wave fieldwave
	
        
	variable numPoints = numpnts(fieldwave)
	variable count = 0
	variable ii = 0

	do
		if (ii >= numPoints)  
			break
		endif
		
		if (fieldwave[ii] != 0)
			break
 		endif
        
		count += 1
		ii += 1
	while (1)
	DeletePoints 0,count, $(kValue + "_time_stamp_rawdata"),$(kValue + "_Temp_chA_rawdata"),$(kValue + "_Temp_chB_rawdata");
	DeletePoints 0,count, $(kValue + "_field_rawdata"),$(kValue + "_voltage_rawdata"),$(kValue + "_resistance_rawdata")
    
end




Function CalcDialog(V_mode)
	variable V_mode
	String CurDataFolder = GetDataFolder(1)
	
	string promptStr = "計算は、 ターゲットフォルダ \""+ CurDataFolder +  "\" 内のすべての温度フォルダに適用されます。\n\n"
	promptStr += "ターゲットフォルダは上記で間違いありませんか？"
	
	Switch (V_mode)
        Case 0:
            promptStr += ""
        break
        
        Case 1:
            promptStr += "\nまた、試料の厚みt (mm)は入力しましたか？"
        break
        
        Case 2:
            promptStr += "\nまた、試料の厚みt (mm)、幅w (mm)、電圧端子間長さL (mm)は入力しましたか？"
        break
    EndSwitch
	
	
	DoAlert 1, promptStr
	
	if (V_Flag != 1)
		Print "Calculation cancelled."
		return -1 // User canceled
	endif
	
	Switch (V_mode)
        Case 0:            
			calc_Ryx_kValues()
        break        
        Case 1:
            calc_rho_yx_kValues()
        break        
        Case 2:
            
        break
    EndSwitch
	
End


Function UpdateWaveList(dfName)
	String dfName
	Variable ii
	
	SVAR wavelists = root:wavelists
	
	String currentDF = GetDataFolder(1)
	if( !DataFolderExists(dfName) )
		Abort dfName+"not found."
	return -1
	endif
	
	SetDataFolder $dfName
	
	String waveNames = WaveList("*K_*",";","")
	
	if (strlen(waveNames) > 0)
		
		Variable numWaves = ItemsInList(waveNames) // Count "*K_*" wave number 
		
		for(ii = 0; ii < numWaves; ii += 1)
			String WaveStr = StringFromList(ii, waveNames)
			String WaveStrChip
			
			variable kIndex = strsearch(WaveStr, "K_", 0)
				
			if (kIndex >= 0)
				WaveStrChip = WaveStr[kIndex + 2, strlen(WaveStr)-1]
			endif
			
			if (strsearch(wavelists, ";" + WaveStrChip + ";", 0) < 0)
				wavelists += WaveStrChip + ";"
			else
			endif
		endfor
	endif

	String folderList = DataFolderDir(1) 
	string processedList = ""
	
	variable startIndex = strsearch(folderList, "FOLDERS:",0) + strlen("FOLDERS:")
	variable endIndex = strsearch(folderList, ";", startIndex)
	
	
	if (startIndex >= 0 && endIndex > startIndex)
		processedList = folderList[startIndex, endIndex - 1]
		
		variable n = 0
		variable folderCount = ItemsInList(processedList, ",")

		do
			if(folderCount == 0)
				break
			endif
			
			string subFolder = StringFromList(n, processedList, ",")
			
			
			if( !DataFolderExists(subFolder) )
			Abort subFolder+"not found."
			return -1
			endif
			
			//Run recursively in all child folders
			setdataFolder $subFolder
			UpdateWaveList(currentDF + "'"+subFolder+"'")
			
			setdataFolder currentDF
			
			n += 1
		while (n < folderCount)
	else
	endif


	if( !DataFolderExists(currentDF) )
		Abort currentDF+"not found."
	return -1
	endif
			
	SetDataFolder currentDF

End




Function PlotDialog()
	String CurDataFolder = GetDataFolder(1)
	
	
	string promptStr = "プロットは、 ターゲットフォルダ \""+ CurDataFolder + "\" 内のすべての温度フォルダに適用されます。\n\n"
	promptStr += "ターゲットフォルダは上記で間違いありませんか？\n\n"
	promptStr += "*注*\ny軸がホール抵抗(率)の場合、測定配置の指定をしてください。配置を間違えた場合は、逆符号のオフセットがかかります。\n"
	promptStr += "また、両軸のwaveの点数が異なる場合は、プロットをスキップします。\n"
	DoAlert 1, promptStr
	
	if (V_Flag != 1)
		Print "Calculation cancelled."
		return -1 // User canceled
	endif
	
	PlotAllWave()
	
End

Function PlotAllWave()
	controlinfo Plot_yVar; String yVar = S_value
	controlinfo Plot_xVar; String xVar = S_value
	
	wave kValues
	
	String currentDF = GetDataFolder(1)
	
	string lstr=""
		
    if (!WaveExists(kValues))
		Print "Error: Expected kValues wave, but it does not exist."
	return 0
	endif

	variable jj = 0
	
	variable kvalue2 = 0
	string kvalue
	
	
		
	Variable YwaveCount = 0
	Variable XwaveCount = 0
		
	do	
		kvalue2 = kvalues[jj]
		
		kvalue = num2str(kvalue2)+"K"
	
		if (!DataFolderExists('kvalue'))
			Print "Error: \"" + kvalue + "\" folder not found."
		return 0
		endif
		
		setdataFolder $kvalue
				
		String YwaveStr = WaveList("*K_"+yVar,";","")
		String XwaveStr = WaveList("*K_"+xVar,";","")
		
		YwaveStr = YwaveStr[0, strlen(YwaveStr)-2]
		XwaveStr = XwaveStr[0, strlen(XwaveStr)-2]
		
		Wave yw = $YwaveStr
		Wave xw = $XwaveStr
        
		if (WaveExists(yw))
			YwaveCount += 1
		endif
        
		if (WaveExists(xw))
			XwaveCount += 1
		endif
		
		if (WaveExists(yw) && WaveExists(xw))
		
      		if (numpnts(yw) != numpnts(xw))
        		SetDataFolder currentDF
				Print "Error: The number of points differs between Ywave and Xwave : '*K_"+yVar+"' vs. '*K_"+xVar+"'"
				Abort "Error: \nThe number of points differs between Ywave and Xwave : '*K_"+yVar+"' vs. '*K_"+xVar+"'"
			endif

			if (jj == 0)
				Display $YwaveStr vs $XwaveStr
				
			
				String Ylabel = GetLabel(YwaveStr)
				String Xlabel = GetLabel(XwaveStr)
				
			else
				AppendToGraph $YwaveStr vs $XwaveStr
			endif
			
			lstr = lstr + "\\s('"+YwaveStr+"') "+num2str(kvalue2)+" K\n"
		
		endif
		
		setdataFolder ::
		jj += 1
	while(jj < numpnts(kValues))
	
	
	if (YwaveCount == 0)
		Print "Error: \"*K_"+ yVar +"\" Wave not found."
		SetDataFolder currentDF
		return -1
	endif
        
	if (XwaveCount == 0)
		Print "Error: \"*K_"+ xVar +"\" Wave not found."
		SetDataFolder currentDF
		return -1
	endif
	
	lstr=lstr[0, strlen(lstr)-2]
	TextBox/B=1/C/N=text0/F=0/A=MC/X=76.5/Y=-7 lstr
	
	Plot_style(Ylabel,Xlabel,YwaveStr)
	
	Print "'*K_"+yVar+"' vs. '*K_"+xVar+"' Plot completed."

End

Function/S GetLabel(axisVar)
	String axisVar
	
	String Output
	
	if (strsearch(axisVar, "square", 0, 2) > 0)
		Output = "\\f02B\\f00\\S2\\M (T\\S2\\M)"
		return Output
	
	elseif (strsearch(axisVar, "field", 0, 2) > 0)
		Output = "\\f02B\\f00 (T)"
		return Output
	
	elseif (strsearch(axisVar, "rho_yx", 0, 2) > 0)
		Output = "\\f02ρ\\Byx\\M\\f00 (Ω cm)"
		return Output
	
	elseif (strsearch(axisVar, "rho_xx", 0, 2) > 0)
		Output = "\\f02ρ\\Bxx\\M\\f00 (Ω cm)"
		return Output
				
	elseif (strsearch(axisVar, "resistance", 0, 2) > 0)
		Output = "Resistance (Ω)"
		return Output
	
	elseif (strsearch(axisVar, "Rxx", 0, 2) > 0)
		Output = "\\f02R\\Bxx\\M\\f00 (Ω) "
		return Output
	
	elseif (strsearch(axisVar, "Ryx", 0, 2) > 0)
		Output = "\\f02R\\Byx\\M\\f00 (Ω) "
		return Output
	
	elseif (strsearch(axisVar, "Temp", 0, 2) > 0)
		Output = "\\f02T\\f00 (K)"
		return Output
	
	elseif (strsearch(axisVar, "MR", 0, 2) > 0)
		Output = "MR (%)"
		return Output
	
	elseif (strsearch(axisVar, "Volt", 0, 2) > 0)
		Output = "Voltage (V)"
		return Output
	else
		return "" // other case
	
	endif
	
	
End



function Plot_style(Ylabel,Xlabel,YwaveStr)
	String Ylabel,Xlabel,YwaveStr
	ModifyGraph/Z margin(right)=141,gfSize=16,width={Aspect,1.2}
	ModifyGraph/Z mode=4
	ModifyGraph/Z marker=19
	ModifyGraph/Z lSize=1.5
	ModifyGraph height=226.772
	ModifyGraph/Z msize=2
	ModifyGraph/Z mrkThick=1.5
	ModifyGraph/Z tick=2
	ModifyGraph/Z zero=4
	ModifyGraph/Z mirror=1
	ModifyGraph/Z standoff=0
	ModifyGraph/Z btLen=7
	ModifyGraph/Z stLen=3
	Label/Z left Ylabel
	Label/Z bottom Xlabel
	
	NVAR gRadioVal= root:gRadioVal
	
	
		
	if ((strsearch(YwaveStr, "Ryx",0,2)>0) || (strsearch(YwaveStr, "rho_yx",0,2)>0) )
	
		if (gRadioVal == 1)
			ModifyGraph muloffset={-1,0}
		elseif (gRadioVal == 2)
			ModifyGraph muloffset={0,0}
		endif
	
	endif
	
	if ((strsearch(YwaveStr, "MR",0,2)>0))
		ModifyGraph muloffset={0,100}
	endif
	
	GraphColorMarker(1, "Turbo",20)
end


function GraphColorMarker(direction, colortab,mVar)
variable direction //up 0, down 1
string colortab
variable mVar

string wl = TraceNameList("",";",1)
string wn
variable ii=0
do
	wn = StringFromList(ii, wl)
	if(strlen(wn) == 0)
		break
	endif
	ii+=1
while(1)
variable mx=ii
ColorTab2Wave $colortab //selest color table turbo may be better
wave M_colors
variable step = floor((numpnts(M_colors))/3/ii)	
ii=0
do
	if(direction==0)
	ModifyGraph rgb($(StringFromList(ii, wl)))=(M_colors[ii*step][0],M_colors[ii*step][1],M_colors[ii*step][2])
	else
	ModifyGraph rgb($(StringFromList(ii, wl)))=(M_colors[(mx-ii-1)*step][0],M_colors[(mx-ii-1)*step][1],M_colors[(mx-ii-1)*step][2])
	endif
	ii+=1
while(ii<mx)
killwaves M_colors

ModifyGraph marker=mVar-1
end

graphMarkers(ColorNum, ColorTab)

/// panel ///

Menu "Magnetic_Field"
"Open Analysis Panel...", ShowMagnetAnlysisPanel()
end

proc ShowMagnetAnlysisPanel()
	DoWindow/F MagnetAnlysisPanel
	if (V_Flag == 0)
		MagnetAnlysisPanel()
	endif
End

Window MagnetAnlysisPanel() : Panel
	Variable/G gRadioVal=root:gRadioVal
	String/G wavelists = root:wavelists
	
	gRadioVal = 1
	wavelists = ";"
	
	PauseUpdate; Silent 1		// building window...	
	NewPanel /N=MagnetAnlysisPanel/K=1 /W=(547,103,1295,830) as "Magnetic Field Measurement Analysis Panel"
	//ShowTools/A
	//ModifyPanel fixedSize=1, noEdit=1
	SetDrawLayer UserBack
	SetDrawEnv linefgc= (48059,48059,48059)
	DrawLine 25,80,403,80
	DrawText 30,174,"\\f03Details\\f00\r   Loads all data that contains the string entered in \\K(65535,0,26214)“Included Text”\\K(0,0,0)\r\n   and \\K(65535,0,26214)“K”\\K(0,0,0) (=Kelvin) in the file name.\r\n\r\n   e.g.) “Included Text” = “Hall”\r\n\t\t\t\t    : 5\\K(65535,0,26214)K\\K(0,0,0)_rampfield_\\K(65535,0,26214)Hall\\K(0,0,0), 10\\K(65535,0,26214)K\\K(0,0,0)_rampfield_\\K(65535,0,26214)Hall\\K(0,0,0), ..."
	DrawText 30,232,"\\f03Steps\\f00\r\n   1. Enter “Include Text” and press “Load \\f02B\\f00 Dependence Data” button.\r\n   2. Select the folder that contains the data you want to load."
	DrawText 30,263,"\\f03*Option*\\f00"
	SetDrawEnv gstart
	SetDrawEnv fillfgc= (56797,56797,56797)
	DrawPoly 81,430,1,1,{275,399,253,414,337,414,359,399,275,399}
	SetDrawEnv fillfgc= (48059,48059,48059)
	DrawPoly 59,473,1,1,{253,442,253,414,337,414,337,442,253,442}
	SetDrawEnv fillfgc= (34952,34952,34952)
	DrawPoly 166,458,1,1,{360,427,337,442,337,414,359,399,360,427}
	SetDrawEnv gstop
	SetDrawEnv linethick= 2,linecap= 1
	DrawLine 131,439,131,407
	SetDrawEnv linethick= 2,linecap= 1
	DrawLine 88,439,88,407
	SetDrawEnv linethick= 2,linecap= 1
	DrawLine 131,407,88,407
	SetDrawEnv linethick= 2
	DrawOval 94,390,124,420
	SetDrawEnv fsize= 14
	DrawText 103,414,"\\f03V\\f00"
	SetDrawEnv fillfgc= (65535,65535,65535,262)
	DrawBezier 58,473,-0.2,1,{237,543,259,543,261,516,237,516}
	SetDrawEnv fillfgc= (65535,65535,65535,262)
	DrawBezier 145.7,473.7,0.4,0.9,{217,536.1,272.1,544.7,296.3,531.4,264,517}
	SetDrawEnv fillfgc= (65535,65535,65535,262)
	DrawBezier 130,441,-0.6,1.5,{211,524,211,532.6,284.8,533.3,284.7,524.4}
	SetDrawEnv linefgc= (48059,48059,48059)
	DrawLine 35,580,391,580
	SetDrawEnv linefgc= (48059,48059,48059)
	DrawLine 725,63,439,63
	SetDrawEnv gstart
	SetDrawEnv gstart
	DrawLine 464.375,142.7,551.375,142.7
	DrawLine 508,120,508,168
	SetDrawEnv linecap= 1,fillfgc= (48059,48059,48059)
	DrawRect 476,127,540,159
	SetDrawEnv fsize= 16,textrgb= (65535,0,26214),textxjust= 1,textyjust= 1
	DrawText 452,141,"\\f02I\\f00+"
	SetDrawEnv fsize= 16,textxjust= 1,textyjust= 1
	DrawText 562,141,"\\f02I\\f00-"
	SetDrawEnv fsize= 16,textxjust= 1,textyjust= 1
	DrawText 507,109,"\\f02V\\f00-"
	SetDrawEnv fsize= 16,textrgb= (65535,0,26214),textxjust= 1,textyjust= 1
	DrawText 507,176,"\\f02V\\f00+"
	SetDrawEnv gstop
	SetDrawEnv gstart
	SetDrawEnv gstart
	SetDrawEnv linethick= 2,fillfgc= (65535,65535,65535,0)
	DrawOval 459,165,474,180
	SetDrawEnv linethick= 2
	DrawLine 461,167,473,179
	SetDrawEnv linethick= 2
	DrawLine 473,167,461,179
	SetDrawEnv gstop
	SetDrawEnv fsize= 16,textxjust= 1,textyjust= 1
	DrawText 449,172,"\\f02B\\f00"
	SetDrawEnv gstop
	SetDrawEnv gstop
	SetDrawEnv gstart
	SetDrawEnv gstart
	DrawLine 611,143,698,143
	DrawLine 655,120,655,167
	SetDrawEnv linecap= 1,fillfgc= (48059,48059,48059)
	DrawRect 623,127,687,159
	SetDrawEnv fsize= 16,textrgb= (65535,0,26214),textxjust= 1,textyjust= 1
	DrawText 599,141,"\\f02I\\f00+"
	SetDrawEnv fsize= 16,textxjust= 1,textyjust= 1
	DrawText 709,141,"\\f02I\\f00-"
	SetDrawEnv fsize= 16,textrgb= (65535,0,26214),textxjust= 1,textyjust= 1
	DrawText 654,109,"\\f02V\\f00+"
	SetDrawEnv fsize= 16,textxjust= 1,textyjust= 1
	DrawText 654,176,"\\f02V\\f00-"
	SetDrawEnv gstop
	SetDrawEnv gstart
	SetDrawEnv gstart
	SetDrawEnv linethick= 2,fillfgc= (65535,65535,65535,0)
	DrawOval 606,165,621,180
	SetDrawEnv linethick= 2
	DrawLine 608,167,620,179
	SetDrawEnv linethick= 2
	DrawLine 620,167,608,179
	SetDrawEnv gstop
	SetDrawEnv fsize= 16,textxjust= 1,textyjust= 1
	DrawText 596,172,"\\f02B\\f00"
	SetDrawEnv gstop
	SetDrawEnv gstop
	DrawText 445,266,"\\f03*Cotion*\\f00\r\n   If the wrong placement is selected, the sign of\r   the Hall resistance will be plotted in reverse."
	SetDrawEnv linefgc= (48059,48059,48059)
	DrawLine 725,323,439,323
	SetDrawEnv fsize= 14,textyjust= 1
	DrawText 573,422,"vs."
	DrawText 445,515,"\\f03*Cotion*\\f00\r\n   If the number of points for\r   Ywave and Xwave are different,\r   the operation will be skipped."
	SetDrawEnv linefgc= (48059,48059,48059)
	DrawLine 725,559,439,559
	SetDrawEnv linefgc= (48059,48059,48059)
	DrawLine 620,602,446,602
	DrawText 445,704,"\\f03Detail\\f00\r\n   Change the plot style of the Active Graph."
	GroupBox LoadGroup,pos={14.00,16.00},size={400.00,306.00},title="Load Files"
	GroupBox LoadGroup,font="Arial",fSize=13,fStyle=0
	SetVariable IncludeText,pos={30.00,47.00},size={189.00,19.00}
	SetVariable IncludeText,title="\\f00Included Text:\\f00"
	SetVariable IncludeText,help={"Enter the name of the data contained in the file name (e.g. Hall, MR,...)"}
	SetVariable IncludeText,font="Arial",fSize=14,fStyle=1
	SetVariable IncludeText,valueColor=(65535,0,26214),value=_STR:"Hall"
	Button loadbtn,pos={231.00,41.00},size={166.00,30.00},proc=ButtonProc
	Button loadbtn,title="Load \\f03B\\f01 Dependence Data..."
	Button loadbtn,help={"Open Folder Selection Dialog"},font="Arial",fStyle=1
	Button loadbtn,fColor=(13107,0,0),valueColor=(65535,65535,65535)
	CheckBox CutValue,pos={38.00,269.00},size={202.00,15.00}
	CheckBox CutValue,title="Cut the first excess zero field data.",font="Arial"
	CheckBox CutValue,value=1
	CheckBox CalcRyxValue,pos={38.00,290.00},size={358.00,17.00}
	CheckBox CalcRyxValue,title=" (Hall measurement) Calc even & odd components (\\f02R\\Bxx\\M\\f00 & \\f02R\\Byx\\M\\f00)."
	CheckBox CalcRyxValue,font="Arial",value=0
	GroupBox AnalysisGroup,pos={14.00,323.00},size={400.00,395.00},title="Analysis"
	GroupBox AnalysisGroup,font="Arial",fSize=13
	TabControl AnalysisTab,pos={19.00,344.00},size={390.00,367.00},proc=TabProc
	TabControl AnalysisTab,font="Arial",fStyle=1,tabLabel(0)="Hall",tabLabel(1)="MR"
	TabControl AnalysisTab,value=0
	GroupBox SampleSizegroup,pos={199.00,376.00},size={185.00,107.00}
	GroupBox SampleSizegroup,title="Sample Size",font="Arial",fSize=12
	SetVariable ThkVar_Alltab,pos={207.00,400.00},size={160.00,18.00},bodyWidth=60
	SetVariable ThkVar_Alltab,title="Thickness (mm) :",font="Arial",format="%.3f"
	SetVariable ThkVar_Alltab,limits={0,100,1},value=_NUM:0
	SetVariable WidVar_DAtab1,pos={233.00,451.00},size={134.00,18.00},bodyWidth=60
	SetVariable WidVar_DAtab1,title="Width (mm) :",font="Arial",format="%.3f"
	SetVariable WidVar_DAtab1,limits={0,100,1},value=_NUM:0
	SetVariable LenVar_DAtab1,pos={226.00,426.00},size={141.00,18.00},bodyWidth=60
	SetVariable LenVar_DAtab1,title="Length (mm) :",font="Arial",format="%.3f"
	SetVariable LenVar_DAtab1,limits={0,100,1},value=_NUM:0
	TitleBox text3,pos={101.00,448.00},size={17.00,19.00},title=" \\f02L\\f00 "
	TitleBox text3,labelBack=(48059,48059,48059),font="Arial",fSize=16,frame=0
	TitleBox text2,pos={162.00,466.00},size={19.00,19.00},title=" \\f02w\\f00 "
	TitleBox text2,labelBack=(63993,63993,63993),font="Arial",fSize=16,frame=0
	TitleBox text1,pos={43.00,451.00},size={12.00,19.00},title=" \\f02t\\f00 "
	TitleBox text1,labelBack=(63993,63993,63993),font="Arial",fSize=16,frame=0
	Button Calc_rho_xx_HDtab1,pos={201.00,493.00},size={182.00,30.00},proc=ButtonProc_7
	Button Calc_rho_xx_HDtab1,title="Calc \\f02ρ\\Bxx\\M\\f00",font="Arial"
	Button Calc_rho_xx_HDtab1,fStyle=0,fColor=(65535,65535,65535)
	Button Calc_MR_HDtab1,pos={201.00,530.00},size={182.00,30.00},proc=ButtonProc_8
	Button Calc_MR_HDtab1,title="Calc MR (%)",font="Arial",fStyle=0
	Button Calc_MR_HDtab1,fColor=(65535,65535,65535)
	Button Calc_Ryx_HDtab0,pos={201.00,493.00},size={182.00,30.00},disable=1,proc=ButtonProc_2
	Button Calc_Ryx_HDtab0,title="Calc \\f02R\\Bxx\\M\\f00 & \\f02R\\Byx\\M\\f00"
	Button Calc_Ryx_HDtab0,font="Arial",fStyle=0,fColor=(65535,65535,65535)
	Button Calc_rho_yx_HDtab0,pos={201.00,530.00},size={182.00,30.00},disable=1,proc=ButtonProc_3
	Button Calc_rho_yx_HDtab0,title="Calc \\f02ρ\\Byx",font="Arial",fStyle=0
	Button Calc_rho_yx_HDtab0,fColor=(65535,65535,65535)
	TitleBox MR_txt1_HDtab1,pos={35.00,591.00},size={344.00,47.00}
	TitleBox MR_txt1_HDtab1,title="\\f03Details\\f00\r\n   You can calculate \\f02ρ\\Bxx\\M\\f00 & MR (%) from the loaded raw data with \r\n   pressing the buttons."
	TitleBox MR_txt1_HDtab1,font="Arial",frame=0,anchor=LC
	TitleBox MR_txt2_HDtab1,pos={35.00,651.00},size={354.00,45.00}
	TitleBox MR_txt2_HDtab1,title="\\f03Steps\\f00\r\n   1. Set the current data folder to the \\K(65535,0,26214)Target Data Folder (e.g. MR)\\K(0,0,0).\r\n   2.  Enter the \\K(65535,0,26214)Thickness (mm), Length (mm), Width (mm)\\K(0,0,0)."
	TitleBox MR_txt2_HDtab1,font="Arial",frame=0,anchor=LC
	TitleBox Hall_txt1_HDtab0,pos={35.00,591.00},size={324.00,47.00},disable=1
	TitleBox Hall_txt1_HDtab0,title="\\f03Details\\f00\r\n   You can calculate \\f02R\\Bxx\\M\\f00, \\f02R\\Byx\\M\\f00, & \\f02ρ\\Byx\\M\\f00 from the loaded raw data \r\n   with pressing the buttons."
	TitleBox Hall_txt1_HDtab0,font="Arial",frame=0,anchor=LC
	TitleBox Hall_txt2_HDtab0,pos={35.00,651.00},size={358.00,47.00},disable=1
	TitleBox Hall_txt2_HDtab0,title="\\f03Steps\\f00\r\n   1. Set the current data folder to the \\K(65535,0,26214)Target Data Folder (e.g. Hall)\\K(0,0,0).\r\n   2. (\\f02ρ\\Byx\\M\\f00 Calc.) Enter the \\K(65535,0,26214)Thickness (mm)\\K(0,0,0) of the sample."
	TitleBox Hall_txt2_HDtab0,font="Arial",frame=0,anchor=LC
	GroupBox group4,pos={429.00,17.00},size={308.00,701.00},title="Plot"
	GroupBox group4,font="Arial",fSize=12
	TitleBox title3,pos={441.00,46.00},size={114.00,15.00}
	TitleBox title3,title="Select Arrangement",font="Arial",frame=0,fStyle=1
	GroupBox RadioTitleBox1,pos={440.00,88.00},size={135.00,104.00},font="Arial"
	TitleBox RadioTitle1,pos={462.00,82.00},size={78.00,15.00},title="Now Selected"
	TitleBox RadioTitle1,labelBack=(61680,61680,61680),font="Arial",fSize=12,frame=0
	TitleBox RadioTitle1,fStyle=1,fColor=(65535,0,26214),anchor=MC
	CheckBox check1,pos={501.00,197.00},size={14.00,14.00},proc=RadioBtnProc
	CheckBox check1,title="",font="Arial",fStyle=1,fColor=(65535,0,26214)
	CheckBox check1,value=1,mode=1
	GroupBox RadioTitleBox2,pos={587.00,88.00},size={135.00,104.00},disable=1
	GroupBox RadioTitleBox2,font="Arial"
	TitleBox RadioTitle2,pos={609.00,82.00},size={78.00,15.00},disable=1
	TitleBox RadioTitle2,title="Now Selected",labelBack=(61680,61680,61680)
	TitleBox RadioTitle2,font="Arial",fSize=12,frame=0,fStyle=1
	TitleBox RadioTitle2,fColor=(65535,0,26214),anchor=MC
	CheckBox check2,pos={648.00,197.00},size={14.00,14.00},proc=RadioBtnProc
	CheckBox check2,title="",font="Arial",fStyle=1,fColor=(65535,0,26214)
	CheckBox check2,value=0,mode=1
	TitleBox title4,pos={441.00,306.00},size={64.00,15.00},title="Plot Waves"
	TitleBox title4,font="Arial",frame=0,fStyle=1
	Button list_update,pos={447.00,345.00},size={270.00,33.00},proc=ButtonProc_5
	Button list_update,title="WaveList  Update",font="Arial",fSize=12,fStyle=0
	Button list_update,fColor=(65535,65535,65535)
	TitleBox title1,pos={447.00,394.00},size={41.00,15.00},title="Ywave :"
	TitleBox title1,font="Arial",fSize=12,frame=0
	TitleBox title2,pos={597.00,394.00},size={41.00,15.00},title="Xwave :"
	TitleBox title2,font="Arial",fSize=12,frame=0
	PopupMenu Plot_yVar,pos={447.00,412.00},size={120.00,19.00},bodyWidth=120,proc=PopMenuProc
	PopupMenu Plot_yVar,font="Arial"
	PopupMenu Plot_yVar,mode=13,popvalue="Ywave",value=#""
	PopupMenu Plot_xVar,pos={597.00,412.00},size={120.00,19.00},bodyWidth=120,proc=PopMenuProc
	PopupMenu Plot_xVar,font="Arial"
	PopupMenu Plot_xVar,mode=9,popvalue="Xwave",value=#""
	Button Plot,pos={637.00,458.00},size={80.00,60.00},proc=ButtonProc_6
	Button Plot,title="Plot !",font="Arial",fSize=14,fStyle=1
	Button Plot,valueColor=(65535,0,26214)
	TitleBox title5,pos={441.00,542.00},size={74.00,15.00},title="Modify Graph"
	TitleBox title5,font="Arial",frame=0,fStyle=1
	PopupMenu MarkerStr,pos={446.00,573.00},size={174.00,19.00},bodyWidth=51
	PopupMenu MarkerStr,title="Set All Markers To:       ",font="Arial"
	PopupMenu MarkerStr,mode=20,value=#"\"*MARKERPOP*\""
	PopupMenu ColorStr,pos={446.00,615.00},size={174.00,19.00},bodyWidth=94
	PopupMenu ColorStr,title="Set Colors To:",font="Arial"
	PopupMenu ColorStr,mode=61,value=#"\"*COLORTABLEPOPNONAMES*\""
	TitleBox title0,pos={447.00,642.00},size={107.00,15.00}
	TitleBox title0,title="Up (0) or Down (1) :",font="Arial",frame=0
	SetVariable ColorVar,pos={564.00,641.00},size={56.00,19.00},title=" "
	SetVariable ColorVar,font="Arial",fSize=14,limits={0,1,1},value=_NUM:1
	Button loadbtn1,pos={637.00,587.00},size={80.00,60.00},proc=ButtonProc_1
	Button loadbtn1,title="Update !",font="Arial",fSize=14,fStyle=1
	Button loadbtn1,valueColor=(65535,0,26214)
	CheckBox Remove_yx_HDtab1,pos={34.00,530.00},size={14.00,14.00},title="",value=0,disable=1
	TitleBox title6_HDtab1,pos={51.00,530.00},size={131.00,30.00},disable=1
	TitleBox title6_HDtab1,title="Remove yx components\nfor MR calculation."
	TitleBox title6_HDtab1,font="Arial",frame=0
EndMacro



Function RadioBtnProc(name,value) : CheckBoxControl
	String name	
	Variable value
	NVAR gRadioVal= root:gRadioVal
	strswitch (name)
	
	case "check1":
		gRadioVal= 1
	break
	case "check2":
		gRadioVal= 2
	break	
	endswitch
		
	if (gRadioVal==1)
		ModifyControl RadioTitle1 disable=0
		ModifyControl RadioTitleBox1 disable=0
		ModifyControl RadioTitle2 disable=1
		ModifyControl RadioTitleBox2 disable=1
		
	elseif  (gRadioVal==2)
		ModifyControl RadioTitle1 disable=1
		ModifyControl RadioTitleBox1 disable=1
		ModifyControl RadioTitle2 disable=0
		ModifyControl RadioTitleBox2 disable=0
	endif
	
	CheckBox check1,value= gRadioVal==1
	CheckBox check2,value= gRadioVal==2
End


Function TabProc(ctrlName,tabNum) : TabControl
	String ctrlName
	Variable tabNum
	
	String DAblTabMatch= "*_DAtab"+num2istr(tabNum)
	String HdTabMatch= "*_HDtab"+num2istr(tabNum)
	String controls= ControlNameList("")
	Variable ii, n= ItemsInList(controls)
	
	for(ii=0; ii<n; ii+=1)
		String control= StringFromList(ii, controls)
		Variable DAblTab= stringmatch(control,"*_DAtab*")
		Variable HdTab= stringmatch(control,"*_HDtab*")
		
		if( DAblTab )
			Variable active= stringmatch(control,DAblTabMatch)
			ModifyControl $control disable=(active ? 0 : 2)
		endif
		if( HdTab )
			Variable show= stringmatch(control,HdTabMatch)
			ModifyControl $control disable=(show ? 0 : 1)
		endif
		
	endfor
	
	

	return 0
End

Function PopMenuProc(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa
	
	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			
			ControlInfo IncludeText
			String Data_value = S_Value
			Load_rampfield_data(data_value)
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function ButtonProc_1(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2:
		
			ControlInfo ColorVar			
			variable ColorNum = V_Value
			
			ControlInfo ColorStr			
			string ColorTab = S_Value
			
			
			ControlInfo MarkerStr
			variable MarkerNum = V_Value
			
			
			GraphColorMarker(ColorNum, ColorTab, MarkerNum)
			
			
			break
		case -1:
			break
	endswitch

	return 0
End

Function ButtonProc_2(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2:
			
			CalcDialog(0)
				
			break
		case -1:
			break
	endswitch

	return 0
End

Function ButtonProc_3(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2:
			
			CalcDialog(1)
				
			break
		case -1:
			break
	endswitch

	return 0
End



Function ButtonProc_5(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2:
			
			SVAR wavelists = root:wavelists
			
			String currentDF = GetDataFolder(1)
			
			SetDataFolder root:			
			wavelists = ";";
			UpdateWaveList("root:")
			
			wavelists = wavelists[1, strlen(wavelists)-2]
			
			SetDataFolder currentDF
			
			String localwavelists = wavelists
			
			localwavelists = "\"" + wavelists + "\""
				
			PopupMenu Plot_yVar value=#localwavelists, mode=1
			PopupMenu Plot_xVar value=#localwavelists, mode=2
			
			ControlUpdate Plot_yVar
			ControlUpdate Plot_xVar
			
			Print "WaveList Update completed."
			break
		case -1:
			break
	endswitch

	return 0
End


Function ButtonProc_6(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2:
		
			PlotDialog()
				
			break
		case -1:
			break
	endswitch

	return 0
End

Function ButtonProc_7(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2:
			
			CalcDialog(2)
			calc_rho_xx_kValues()
				
			break
		case -1:
			break
	endswitch

	return 0
End

Function ButtonProc_8(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2:
			
			CalcDialog(2)
			calc_MR_kValues()
				
			break
		case -1:
			break
	endswitch

	return 0
End
