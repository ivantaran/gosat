package gosat

type elsetrec struct {
	satnum        int // TODO: long int originally
	epochyr       int
	epochtynumrev int
	error         int
	operationmode rune
	init          rune
	method        rune

	/* Near Earth */
	isimp   int
	aycof   float64
	con41   float64
	cc1     float64
	cc4     float64
	cc5     float64
	d2      float64
	d3      float64
	d4      float64
	delmo   float64
	eta     float64
	argpdot float64
	omgcof  float64
	sinmao  float64
	t       float64
	t2cof   float64
	t3cof   float64
	t4cof   float64
	t5cof   float64
	x1mth2  float64
	x7thm1  float64
	mdot    float64
	nodedot float64
	xlcof   float64
	xmcof   float64
	nodecf  float64

	/* Deep Space */
	irez  int
	d2201 float64
	d2211 float64
	d3210 float64
	d3222 float64
	d4410 float64
	d4422 float64
	d5220 float64
	d5232 float64
	d5421 float64
	d5433 float64
	dedt  float64
	del1  float64
	del2  float64
	del3  float64
	didt  float64
	dmdt  float64
	dnodt float64
	domdt float64
	e3    float64
	ee2   float64
	peo   float64
	pgho  float64
	pho   float64
	pinco float64
	plo   float64
	se2   float64
	se3   float64
	sgh2  float64
	sgh3  float64
	sgh4  float64
	sh2   float64
	sh3   float64
	si2   float64
	si3   float64
	sl2   float64
	sl3   float64
	sl4   float64
	gsto  float64
	xfact float64
	xgh2  float64
	xgh3  float64
	xgh4  float64
	xh2   float64
	xh3   float64
	xi2   float64
	xi3   float64
	xl2   float64
	xl3   float64
	xl4   float64
	xlamo float64
	zmol  float64
	zmos  float64
	atime float64
	xli   float64
	xni   float64
	/*  */
	a           float64
	altp        float64
	alta        float64
	epochdays   float64
	jdsatepoch  float64
	jdsatepochF float64
	nddot       float64
	ndot        float64
	bstar       float64
	rcse        float64
	inclo       float64
	nodeo       float64
	ecco        float64
	argpo       float64
	mo          float64
	no_kozai    float64

	// sgp4fix add new variables from tle
	classification rune
	intldesg       string // TODO: length must be 11
	ephtype        int
	elnum          int // TODO: long originally
	revnum         int // TODO: long originally
	// sgp4fix add unkozai'd variable
	no_unkozai float64
	// sgp4fix add singly averaged variables
	am float64
	em float64
	im float64
	Om float64
	om float64
	mm float64
	nm float64
	// sgp4fix add constant parameters to eliminate mutliple calls during execution
	tumin         float64
	mu            float64
	radiusearthkm float64
	xke           float64
	j2            float64
	j3            float64
	j4            float64
	j3oj2         float64

	// Additional elements to capture relevant TLE and object information:
	dia_mm      int     //RSO dia in mm // TODO: long originally
	period_sec  float64 //Period in seconds
	active      bool    //"Active S/C" flag (0=n, 1=y)
	not_orbital bool    //"Orbiting S/C" flag (0=n, 1=y)
	rcs_m2      float64 //"RCS (m^2)" storage
}
