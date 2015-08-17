

./Bin/create-model-lame.x:   lin7@cees-gpu   Mon Aug  3 15:24:51 2015
		n1=300
		o1=0
		d1=4
		n2=300
		o2=0
		d2=4
		label1="z"
		label2="x"
		sets next: in="/data/lin7/rho-layer.H@"
		data_format="xdr_float"


Smooth:   lin7@cees-gpu   Mon Aug  3 15:24:51 2015
		sets next: in="/data/lin7/Smooth@yD1ZBu"
	Data smoothed with tridiagonal solver leaky().
		Data smoothed with: rect1=20
		Data smoothed with: rect2=20
		Data smoothed with: rect3=1
		data_format="xdr_float"

