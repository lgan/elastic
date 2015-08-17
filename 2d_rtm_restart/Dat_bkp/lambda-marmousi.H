No input header 

Su2sep:   gcalves@gertie   Mon May 11 13:24:28 2015
		sets next: in="/net/gertie/scr1/gcalves/vp_marmousi-ii.H@"
		junkME=0
		hff="vp_marmousi-ii.H@@"
		same_record_number=1
	 	n2=557121536 o2=1 d2=1 label2="trace number" 

n3=1 n4=1 n5=1 n6=1 n7=1 n8=1 n9=1 

	n1=2801  o1=0.000000  d1=0.004816   label1="none"   unit1="none"
		gff="-1"
		esize=4
		data_format="xdr_float"
		junkME=0
		same_record_number=1
	 	n2=13601 o2=1 d2=1.25 label2="x" unit2="m"

n3=1 n4=1 n5=1 n6=1 n7=1 n8=1 n9=1 

	n1=2801  o1=0.000000  d1=1.25   label1="z"   unit1="m"
		gff="-1"
		esize=4
hff=-1


/opt/SEP/SEP7.0/bin/Math:   gcalves@gertie   Mon May 11 17:06:04 2015
		sets next: in="/net/gertie/scr1/gcalves/Math@MYOe5S"
		esize=4
	#Using file3*(file1*file1-2.0*file2*file2)*1000000000.0 to produce real data

		#file1="vp_marmousi-ii.H"
		#file2="vs_marmousi-ii.H"
		#file3="density_marmousi-ii.H"
		data_format="xdr_float"
