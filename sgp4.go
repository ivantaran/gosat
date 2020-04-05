package gosat

import "math"

//Constants
const (
	TwoPi   = math.Pi * 2.0
	Deg2Rad = math.Pi / 180.0
	X2O3    = 2.0 / 3.0
	Temp4   = 1.5e-12
)

type elsetrec struct {
	satnum        int // TODO: long int originally
	epochyr       int
	epochtynumrev int
	error         int
	operationmode rune
	init          rune
	method        rune

	/* Near Earth */
	isimp   int // TODO: is bool?
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

type tle struct {
	satn   int
	epoch  float64
	xbstar float64
	xndot  float64
	xnddot float64
	xecco  float64
	xargpo float64
	xinclo float64
	xmo    float64
	xno    float64
	xnodeo float64
}

/* -----------------------------------------------------------------------------
*
*                           procedure dpper
*
*  this procedure provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep        - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
func (s *elsetrec) dpper() (ep float64, inclp float64, nodep float64, argpp float64, mp float64) {
	/* ---------------------- constants ----------------------------- */
	const (
		zns = 1.19459e-5
		zes = 0.01675
		znl = 1.5835218e-4
		zel = 0.05490
	)
	/* --------------- calculate time varying periodics ----------- */
	zm := s.zmos + zns*s.t
	// be sure that the initial call has time set to zero
	if s.init == 'y' { // TODO: else..
		zm = s.zmos
	}
	zf := zm + 2.0*zes*math.Sin(zm)
	sinzf := math.Sin(zf)
	f2 := 0.5*sinzf*sinzf - 0.25
	f3 := -0.5 * sinzf * math.Cos(zf)
	ses := s.se2*f2 + s.se3*f3
	sis := s.si2*f2 + s.si3*f3
	sls := s.sl2*f2 + s.sl3*f3 + s.sl4*sinzf
	sghs := s.sgh2*f2 + s.sgh3*f3 + s.sgh4*sinzf
	shs := s.sh2*f2 + s.sh3*f3
	zm = s.zmol + znl*s.t
	if s.init == 'y' { // TODO: else..
		zm = s.zmol
	}
	zf = zm + 2.0*zel*math.Sin(zm)
	sinzf = math.Sin(zf)
	f2 = 0.5*sinzf*sinzf - 0.25
	f3 = -0.5 * sinzf * math.Cos(zf)
	sel := s.ee2*f2 + s.e3*f3
	sil := s.xi2*f2 + s.xi3*f3
	sll := s.xl2*f2 + s.xl3*f3 + s.xl4*sinzf
	sghl := s.xgh2*f2 + s.xgh3*f3 + s.xgh4*sinzf
	shll := s.xh2*f2 + s.xh3*f3
	pe := ses + sel
	pinc := sis + sil
	pl := sls + sll
	pgh := sghs + sghl
	ph := shs + shll

	if s.init == 'n' {
		pe = pe - s.peo
		pinc = pinc - s.pinco
		pl = pl - s.plo
		pgh = pgh - s.pgho
		ph = ph - s.pho
		inclp = inclp + pinc
		ep = ep + pe
		sinip := math.Sin(inclp)
		cosip := math.Cos(inclp)

		/* ----------------- apply periodics directly ------------ */
		//  sgp4fix for lyddane choice
		//  strn3 used original inclination - this is technically feasible
		//  gsfc used perturbed inclination - also technically feasible
		//  probably best to readjust the 0.2 limit value and limit discontinuity
		//  0.2 rad = 11.45916 deg
		//  use next line for original strn3 approach and original inclination
		//  if (inclo >= 0.2)
		//  use next line for gsfc version and perturbed inclination
		if inclp >= 0.2 {
			ph = ph / sinip
			pgh = pgh - cosip*ph
			argpp = argpp + pgh
			nodep = nodep + ph
			mp = mp + pl
		} else {
			/* ---- apply periodics with lyddane modification ---- */
			sinop := math.Sin(nodep)
			cosop := math.Cos(nodep)
			alfdp := sinip * sinop
			betdp := sinip * cosop
			dalf := ph*cosop + pinc*cosip*sinop
			dbet := -ph*sinop + pinc*cosip*cosop
			alfdp = alfdp + dalf
			betdp = betdp + dbet
			nodep = math.Mod(nodep, TwoPi)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (nodep < 0.0) && (s.operationmode == 'a') {
				nodep = nodep + TwoPi
			}
			xls := mp + argpp + cosip*nodep
			dls := pl + pgh - pinc*nodep*sinip
			xls = xls + dls
			xnoh := nodep
			nodep = math.Atan2(alfdp, betdp)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (nodep < 0.0) && (s.operationmode == 'a') {
				nodep = nodep + TwoPi
			}
			if math.Abs(xnoh-nodep) > math.Pi {
				if nodep < xnoh {
					nodep = nodep + TwoPi
				} else {
					nodep = nodep - TwoPi
				}
			}
			mp = mp + pl
			argpp = xls - mp - cosip*nodep
		}
	}
	return
}

/*-----------------------------------------------------------------------------
*
*                           procedure initl
*
*  this procedure initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satn        - satellite number - not needed, placed in satrec
*    xke         - reciprocal of tumin
*    j2          - j2 zonal harmonic
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    getgravconst- no longer used
*    gstime      - find greenwich sidereal time from the julian date
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func (s *elsetrec) initl(epoch float64) (
	ainv float64, ao float64, cosio float64, cosio2 float64, con42 float64,
	sinio float64, omeosq float64, posq float64, rp float64, rteosq float64) {

	/* ------------- calculate auxillary epoch quantities ---------- */
	eccsq := s.ecco * s.ecco
	omeosq = 1.0 - eccsq
	rteosq = math.Sqrt(omeosq)
	cosio = math.Cos(s.inclo)
	cosio2 = cosio * cosio

	/* ------------------ un-kozai the mean motion ----------------- */
	ak := math.Pow(s.xke/s.no_kozai, X2O3)
	d1 := 0.75 * s.j2 * (3.0*cosio2 - 1.0) / (rteosq * omeosq)
	del := d1 / (ak * ak)
	adel := ak * (1.0 - del*del - del*(1.0/3.0+134.0*del*del/81.0))
	del = d1 / (adel * adel)
	s.no_unkozai = s.no_kozai / (1.0 + del)

	ao = math.Pow(s.xke/s.no_unkozai, X2O3)
	sinio = math.Sin(s.inclo)
	po := ao * omeosq
	con42 = 1.0 - 5.0*cosio2
	s.con41 = -con42 - cosio2 - cosio2
	ainv = 1.0 / ao
	posq = po * po
	rp = ao * (1.0 - s.ecco)
	s.method = 'n'

	ts70 := epoch - 7305.0
	ds70 := math.Floor(ts70 + 1.0e-8)
	tfrac := ts70 - ds70
	// find greenwich location at epoch
	c1 := 1.72027916940703639e-2
	thgr70 := 1.7321343856509374
	fk5r := 5.07551419432269442e-15
	c1p2p := c1 + TwoPi
	gsto1 := math.Mod(thgr70+c1*ds70+c1p2p*tfrac+ts70*ts70*fk5r, TwoPi)
	if gsto1 < 0.0 {
		gsto1 = gsto1 + TwoPi
	}
	s.gsto = gstime(epoch + 2433281.5)

	return
}

func gstime(jdut1 float64) float64 {

	tut1 := (jdut1 - 2451545.0) / 36525.0
	temp := -6.2e-6*tut1*tut1*tut1 +
		0.093104*tut1*tut1 +
		(876600.0*3600+8640184.812866)*tut1 + 67310.54841 // seconds
	temp = math.Mod(temp*Deg2Rad/240.0, TwoPi) //360/86400 = 1/240, to deg, to rad

	// ------------------------ check quadrants ---------------------
	if temp < 0.0 {
		temp += TwoPi
	}

	return temp
}

/* -----------------------------------------------------------------------------
*
*                           function getgravconst
*
*  this function gets constants for the propagator. note that mu is identified to
*    facilitiate comparisons with newer models. the common useage is wgs72.
*
*  author        : david vallado                  719-573-2600   21 jul 2006
*
*  inputs        :
*    whichconst  - which set of constants to use  wgs72old, wgs72, wgs84
*
*  outputs       :
*    tumin       - minutes in one time unit
*    mu          - earth gravitational parameter
*    radiusearthkm - radius of the earth in km
*    xke         - reciprocal of tumin
*    j2, j3, j4  - un-normalized zonal harmonic values
*    j3oj2       - j3 divided by j2
*
*  locals        :
*
*  coupling      :
*    none
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
--------------------------------------------------------------------------- */

func (s *elsetrec) getgravconst(whichconst string) {
	switch whichconst {
	// -- wgs-72 low precision str#3 constants --
	case "wgs72old":
		s.mu = 398600.79964        // in km3 / s2
		s.radiusearthkm = 6378.135 // km
		s.xke = 0.0743669161       // reciprocal of tumin
		s.tumin = 1.0 / s.xke
		s.j2 = 0.001082616
		s.j3 = -0.00000253881
		s.j4 = -0.00000165597
		s.j3oj2 = s.j3 / s.j2
	// ------------ wgs-72 constants ------------
	case "wgs72":
		s.mu = 398600.8            // in km3 / s2
		s.radiusearthkm = 6378.135 // km
		s.xke = 60.0 / math.Sqrt(s.radiusearthkm*s.radiusearthkm*s.radiusearthkm/s.mu)
		s.tumin = 1.0 / s.xke
		s.j2 = 0.001082616
		s.j3 = -0.00000253881
		s.j4 = -0.00000165597
		s.j3oj2 = s.j3 / s.j2
	case "wgs84":
	default:
		// ------------ wgs-84 constants ------------
		s.mu = 398600.5            // in km3 / s2
		s.radiusearthkm = 6378.137 // km
		s.xke = 60.0 / math.Sqrt(s.radiusearthkm*s.radiusearthkm*s.radiusearthkm/s.mu)
		s.tumin = 1.0 / s.xke
		s.j2 = 0.00108262998905
		s.j3 = -0.00000253215306
		s.j4 = -0.00000161098761
		s.j3oj2 = s.j3 / s.j2
	}
}

/*-----------------------------------------------------------------------------
*
*                           procedure dscom
*
*  this procedure provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
func (s *elsetrec) dscom(
	epoch float64, tc float64,
) (
	snodm float64, cnodm float64, sinim float64, cosim float64, sinomm float64, cosomm float64, day float64,
	em float64, emsq float64, gam float64, rtemsq float64,
	s1 float64, s2 float64, s3 float64, s4 float64, s5 float64, s6 float64, s7 float64,
	ss1 float64, ss2 float64, ss3 float64, ss4 float64, ss5 float64, ss6 float64, ss7 float64,
	sz1 float64, sz2 float64, sz3 float64, sz11 float64, sz12 float64, sz13 float64,
	sz21 float64, sz22 float64, sz23 float64, sz31 float64, sz32 float64, sz33 float64,
	nm float64,
	z1 float64, z2 float64, z3 float64, z11 float64, z12 float64, z13 float64,
	z21 float64, z22 float64, z23 float64, z31 float64, z32 float64, z33 float64,
) {

	const (
		zes    = 0.01675
		zel    = 0.05490
		c1ss   = 2.9864797e-6
		c1l    = 4.7968065e-7
		zsinis = 0.39785416
		zcosis = 0.91744867
		zcosgs = 0.1945905
		zsings = -0.98088458
	)

	nm = s.no_unkozai
	em = s.ecco
	snodm = math.Sin(s.nodeo)
	cnodm = math.Cos(s.nodeo)
	sinomm = math.Sin(s.argpo)
	cosomm = math.Cos(s.argpo)
	sinim = math.Sin(s.inclo)
	cosim = math.Cos(s.inclo)
	emsq = em * em
	betasq := 1.0 - emsq
	rtemsq = math.Sqrt(betasq)

	/* ----------------- initialize lunar solar terms --------------- */
	s.peo = 0.0
	s.pinco = 0.0
	s.plo = 0.0
	s.pgho = 0.0
	s.pho = 0.0
	day = epoch + 18261.5 + tc/1440.0
	xnodce := math.Mod(4.5236020-9.2422029e-4*day, TwoPi)
	stem := math.Sin(xnodce)
	ctem := math.Cos(xnodce)
	zcosil := 0.91375164 - 0.03568096*ctem
	zsinil := math.Sqrt(1.0 - zcosil*zcosil)
	zsinhl := 0.089683511 * stem / zsinil
	zcoshl := math.Sqrt(1.0 - zsinhl*zsinhl)
	gam = 5.8351514 + 0.0019443680*day
	zx := 0.39785416 * stem / zsinil
	zy := zcoshl*ctem + 0.91744867*zsinhl*stem
	zx = math.Atan2(zx, zy)
	zx = gam + zx - xnodce
	zcosgl := math.Cos(zx)
	zsingl := math.Sin(zx)

	/* ------------------------- do solar terms --------------------- */
	zcosg := zcosgs
	zsing := zsings
	zcosi := zcosis
	zsini := zsinis
	zcosh := cnodm
	zsinh := snodm
	cc := c1ss
	xnoi := 1.0 / nm

	for lsflg := 1; lsflg <= 2; lsflg++ {
		a1 := zcosg*zcosh + zsing*zcosi*zsinh
		a3 := -zsing*zcosh + zcosg*zcosi*zsinh
		a7 := -zcosg*zsinh + zsing*zcosi*zcosh
		a8 := zsing * zsini
		a9 := zsing*zsinh + zcosg*zcosi*zcosh
		a10 := zcosg * zsini
		a2 := cosim*a7 + sinim*a8
		a4 := cosim*a9 + sinim*a10
		a5 := -sinim*a7 + cosim*a8
		a6 := -sinim*a9 + cosim*a10

		x1 := a1*cosomm + a2*sinomm
		x2 := a3*cosomm + a4*sinomm
		x3 := -a1*sinomm + a2*cosomm
		x4 := -a3*sinomm + a4*cosomm
		x5 := a5 * sinomm
		x6 := a6 * sinomm
		x7 := a5 * cosomm
		x8 := a6 * cosomm

		z31 = 12.0*x1*x1 - 3.0*x3*x3
		z32 = 24.0*x1*x2 - 6.0*x3*x4
		z33 = 12.0*x2*x2 - 3.0*x4*x4
		z1 = 3.0*(a1*a1+a2*a2) + z31*emsq
		z2 = 6.0*(a1*a3+a2*a4) + z32*emsq
		z3 = 3.0*(a3*a3+a4*a4) + z33*emsq
		z11 = -6.0*a1*a5 + emsq*(-24.0*x1*x7-6.0*x3*x5)
		z12 = -6.0*(a1*a6+a3*a5) + emsq*(-24.0*(x2*x7+x1*x8)-6.0*(x3*x6+x4*x5))
		z13 = -6.0*a3*a6 + emsq*(-24.0*x2*x8-6.0*x4*x6)
		z21 = 6.0*a2*a5 + emsq*(24.0*x1*x5-6.0*x3*x7)
		z22 = 6.0*(a4*a5+a2*a6) + emsq*(24.0*(x2*x5+x1*x6)-6.0*(x4*x7+x3*x8))
		z23 = 6.0*a4*a6 + emsq*(24.0*x2*x6-6.0*x4*x8)
		z1 = z1 + z1 + betasq*z31
		z2 = z2 + z2 + betasq*z32
		z3 = z3 + z3 + betasq*z33
		s3 = cc * xnoi
		s2 = -0.5 * s3 / rtemsq
		s4 = s3 * rtemsq
		s1 = -15.0 * em * s4
		s5 = x1*x3 + x2*x4
		s6 = x2*x3 + x1*x4
		s7 = x2*x4 - x1*x3

		/* ----------------------- do lunar terms ------------------- */
		if lsflg == 1 {
			ss1 = s1
			ss2 = s2
			ss3 = s3
			ss4 = s4
			ss5 = s5
			ss6 = s6
			ss7 = s7
			sz1 = z1
			sz2 = z2
			sz3 = z3
			sz11 = z11
			sz12 = z12
			sz13 = z13
			sz21 = z21
			sz22 = z22
			sz23 = z23
			sz31 = z31
			sz32 = z32
			sz33 = z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl*cnodm + zsinhl*snodm
			zsinh = snodm*zcoshl - cnodm*zsinhl
			cc = c1l
		}
	}

	zmol := math.Mod(4.7199672+0.22997150*day-gam, TwoPi)
	zmos := math.Mod(6.2565837+0.017201977*day, TwoPi)

	/* ------------------------ do solar terms ---------------------- */
	se2 := 2.0 * ss1 * ss6
	se3 := 2.0 * ss1 * ss7
	si2 := 2.0 * ss2 * sz12
	si3 := 2.0 * ss2 * (sz13 - sz11)
	s.sl2 = -2.0 * ss3 * sz2
	s.sl3 = -2.0 * ss3 * (sz3 - sz1)
	s.sl4 = -2.0 * ss3 * (-21.0 - 9.0*emsq) * zes
	s.sgh2 = 2.0 * ss4 * sz32
	s.sgh3 = 2.0 * ss4 * (sz33 - sz31)
	s.sgh4 = -18.0 * ss4 * zes
	s.sh2 = -2.0 * ss2 * sz22
	s.sh3 = -2.0 * ss2 * (sz23 - sz21)

	/* ------------------------ do lunar terms ---------------------- */
	s.ee2 = 2.0 * s1 * s6
	s.e3 = 2.0 * s1 * s7
	s.xi2 = 2.0 * s2 * z12
	s.xi3 = 2.0 * s2 * (z13 - z11)
	s.xl2 = -2.0 * s3 * z2
	s.xl3 = -2.0 * s3 * (z3 - z1)
	s.xl4 = -2.0 * s3 * (-21.0 - 9.0*emsq) * zel
	s.xgh2 = 2.0 * s4 * z32
	s.xgh3 = 2.0 * s4 * (z33 - z31)
	s.xgh4 = -18.0 * s4 * zel
	s.xh2 = -2.0 * s2 * z22
	s.xh3 = -2.0 * s2 * (z23 - z21)

	return
}

func (s *elsetrec) sgp4init(whichconst string, t *tle, opsmode rune, satn int, epoch float64) {
	/* ----------- set all near earth variables to zero ------------ */
	s.isimp = 0
	s.method = 'n'
	s.aycof = 0.0
	s.con41 = 0.0
	s.cc1 = 0.0
	s.cc4 = 0.0
	s.cc5 = 0.0
	s.d2 = 0.0
	s.d3 = 0.0
	s.d4 = 0.0
	s.delmo = 0.0
	s.eta = 0.0
	s.argpdot = 0.0
	s.omgcof = 0.0
	s.sinmao = 0.0
	s.t = 0.0
	s.t2cof = 0.0
	s.t3cof = 0.0
	s.t4cof = 0.0
	s.t5cof = 0.0
	s.x1mth2 = 0.0
	s.x7thm1 = 0.0
	s.mdot = 0.0
	s.nodedot = 0.0
	s.xlcof = 0.0
	s.xmcof = 0.0
	s.nodecf = 0.0

	/* ----------- set all deep space variables to zero ------------ */
	s.irez = 0
	s.d2201 = 0.0
	s.d2211 = 0.0
	s.d3210 = 0.0
	s.d3222 = 0.0
	s.d4410 = 0.0
	s.d4422 = 0.0
	s.d5220 = 0.0
	s.d5232 = 0.0
	s.d5421 = 0.0
	s.d5433 = 0.0
	s.dedt = 0.0
	s.del1 = 0.0
	s.del2 = 0.0
	s.del3 = 0.0
	s.didt = 0.0
	s.dmdt = 0.0
	s.dnodt = 0.0
	s.domdt = 0.0
	s.e3 = 0.0
	s.ee2 = 0.0
	s.peo = 0.0
	s.pgho = 0.0
	s.pho = 0.0
	s.pinco = 0.0
	s.plo = 0.0
	s.se2 = 0.0
	s.se3 = 0.0
	s.sgh2 = 0.0
	s.sgh3 = 0.0
	s.sgh4 = 0.0
	s.sh2 = 0.0
	s.sh3 = 0.0
	s.si2 = 0.0
	s.si3 = 0.0
	s.sl2 = 0.0
	s.sl3 = 0.0
	s.sl4 = 0.0
	s.gsto = 0.0
	s.xfact = 0.0
	s.xgh2 = 0.0
	s.xgh3 = 0.0
	s.xgh4 = 0.0
	s.xh2 = 0.0
	s.xh3 = 0.0
	s.xi2 = 0.0
	s.xi3 = 0.0
	s.xl2 = 0.0
	s.xl3 = 0.0
	s.xl4 = 0.0
	s.xlamo = 0.0
	s.zmol = 0.0
	s.zmos = 0.0
	s.atime = 0.0
	s.xli = 0.0
	s.xni = 0.0

	/* ------------------------ earth constants ----------------------- */
	s.getgravconst(whichconst)

	s.error = 0
	s.operationmode = opsmode
	s.satnum = satn

	s.bstar = t.xbstar
	s.ndot = t.xndot
	s.nddot = t.xnddot
	s.ecco = t.xecco
	s.argpo = t.xargpo
	s.inclo = t.xinclo
	s.mo = t.xmo
	s.no_kozai = t.xno
	s.nodeo = t.xnodeo

	s.am = 0.0
	s.em = 0.0
	s.im = 0.0
	s.Om = 0.0
	s.mm = 0.0
	s.nm = 0.0

	ss := 78.0/s.radiusearthkm + 1.0
	qzms2ttemp := (120.0 - 78.0) / s.radiusearthkm
	qzms2t := qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp

	s.init = 'y'
	s.t = 0.0

	ainv, ao, cosio, cosio2, con42, sinio, omeosq, posq, rp, rteosq := s.initl(epoch)

	s.a = math.Pow(s.no_unkozai*s.tumin, -X2O3)
	s.alta = s.a*(1.0+s.ecco) - 1.0
	s.altp = s.a*(1.0-s.ecco) - 1.0
	s.error = 0

	if (omeosq >= 0.0) || (s.no_unkozai >= 0.0) {
		s.isimp = 0
		if rp < (220.0/s.radiusearthkm + 1.0) {
			s.isimp = 1
		}
		sfour := ss
		qzms24 := qzms2t
		perige := (rp - 1.0) * s.radiusearthkm

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if perige < 156.0 {
			sfour = perige - 78.0
			if perige < 98.0 {
				sfour = 20.0
			}
			qzms24temp := (120.0 - sfour) / s.radiusearthkm
			qzms24 := qzms24temp * qzms24temp * qzms24temp * qzms24temp
			sfour = sfour/s.radiusearthkm + 1.0
		}
		pinvsq := 1.0 / posq

		tsi := 1.0 / (ao - sfour)
		s.eta = ao * s.ecco * tsi
		etasq := s.eta * s.eta
		eeta := s.ecco * s.eta
		psisq := math.Abs(1.0 - etasq)
		coef := qzms24 * math.Pow(tsi, 4.0)
		coef1 := coef / math.Pow(psisq, 3.5)
		cc2 := coef1 * s.no_unkozai * (ao*(1.0+1.5*etasq+eeta*(4.0+etasq)) +
			0.375*s.j2*tsi/psisq*s.con41*(8.0+3.0*etasq*(8.0+etasq)))
		s.cc1 = s.bstar * cc2
		cc3 := 0.0
		if s.ecco > 1.0e-4 {
			cc3 = -2.0 * coef * tsi * s.j3oj2 * s.no_unkozai * sinio / s.ecco
		}
		s.x1mth2 = 1.0 - cosio2
		s.cc4 = 2.0 * s.no_unkozai * coef1 * ao * omeosq * (s.eta*(2.0+0.5*etasq) +
			s.ecco*(0.5+2.0*etasq) -
			s.j2*tsi/(ao*psisq)*(-3.0*s.con41*(1.0-2.0*eeta+etasq*(1.5-0.5*eeta))+
				0.75*s.x1mth2*(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*s.argpo)))
		s.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75*(etasq+eeta) + eeta*etasq)
		cosio4 := cosio2 * cosio2
		temp1 := 1.5 * s.j2 * pinvsq * s.no_unkozai
		temp2 := 0.5 * temp1 * s.j2 * pinvsq
		temp3 := -0.46875 * s.j4 * pinvsq * pinvsq * s.no_unkozai
		s.mdot = s.no_unkozai + 0.5*temp1*rteosq*s.con41 +
			0.0625*temp2*rteosq*(13.0-78.0*cosio2+137.0*cosio4)
		s.argpdot = -0.5*temp1*con42 + 0.0625*temp2*(7.0-114.0*cosio2+395.0*cosio4) +
			temp3*(3.0-36.0*cosio2+49.0*cosio4)
		xhdot1 := -temp1 * cosio
		s.nodedot = xhdot1 + (0.5*temp2*(4.0-19.0*cosio2)+2.0*temp3*(3.0-7.0*cosio2))*cosio
		xpidot := s.argpdot + s.nodedot
		s.omgcof = s.bstar * cc3 * math.Cos(s.argpo)
		s.xmcof = 0.0
		if s.ecco > 1.0e-4 {
			s.xmcof = -X2O3 * coef * s.bstar / eeta
		}
		s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1
		s.t2cof = 1.5 * s.cc1

		if math.Abs(cosio+1.0) > 1.5e-12 {
			s.xlcof = -0.25 * s.j3oj2 * sinio * (3.0 + 5.0*cosio) / (1.0 + cosio)
		} else {
			s.xlcof = -0.25 * s.j3oj2 * sinio * (3.0 + 5.0*cosio) / Temp4
		}
		s.aycof = -0.5 * s.j3oj2 * sinio
		delmotemp := 1.0 + s.eta*math.Cos(s.mo)
		s.delmo = delmotemp * delmotemp * delmotemp
		s.sinmao = math.Sin(s.mo)
		s.x7thm1 = 7.0*cosio2 - 1.0

		/* --------------- deep space initialization ------------- */
		if (2 * math.Pi / s.no_unkozai) >= 225.0 {
			s.method = 'd'
			s.isimp = 1
			tc := 0.0
			inclm := s.inclo

			s.dscom(epoch, tc)

			dpper(s.e3, s.ee2, s.peo, s.pgho, s.pho, s.pinco,
				s.plo, s.se2, s.se3, s.sgh2, s.sgh3, s.sgh4,
				s.sh2, s.sh3, s.si2, s.si3, s.sl2, s.sl3,
				s.sl4, s.t, s.xgh2, s.xgh3, s.xgh4, s.xh2,
				s.xh3, s.xi2, s.xi3, s.xl2, s.xl3, s.xl4,
				s.zmol, s.zmos, inclm, s.init, s.ecco, s.inclo,
				s.nodeo, s.argpo, s.mo, s.operationmode)

			argpm := 0.0
			nodem := 0.0
			mm := 0.0

			dsinit(s.xke, cosim, emsq, s.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3,
				ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, s.t, tc,
				s.gsto, s.mo, s.mdot, s.no_unkozai, s.nodeo,
				s.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33, s.ecco, eccsq,
				em, argpm, inclm, mm, nm, nodem, s.irez, s.atime, s.d2201,
				s.d2211, s.d3210, s.d3222, s.d4410, s.d4422,
				s.d5220, s.d5232, s.d5421, s.d5433, s.dedt, s.didt,
				s.dmdt, dndt, s.dnodt, s.domdt, s.del1, s.del2,
				s.del3, s.xfact, s.xlamo, s.xli, s.xni)
		}

		/* ----------- set variables if not deep space ----------- */
		if s.isimp != 1 {
			cc1sq = s.cc1 * s.cc1
			s.d2 = 4.0 * ao * tsi * cc1sq
			temp = s.d2 * tsi * s.cc1 / 3.0
			s.d3 = (17.0*ao + sfour) * temp
			s.d4 = 0.5 * temp * ao * tsi * (221.0*ao + 31.0*sfour) * s.cc1
			s.t3cof = s.d2 + 2.0*cc1sq
			s.t4cof =
				0.25 * (3.0*s.d3 + s.cc1*(12.0*s.d2+10.0*cc1sq))
			s.t5cof =
				0.2 * (3.0*s.d4 + 12.0*s.cc1*s.d3 +
					6.0*s.d2*s.d2 + 15.0*cc1sq*(2.0*s.d2+cc1sq))
		}
	} // if omeosq = 0 ...

}
