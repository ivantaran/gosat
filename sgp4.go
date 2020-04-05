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
	ainv float64, sinio float64, posq float64, rp float64) {

	/* ------------- calculate auxillary epoch quantities ---------- */
	eccsq := s.ecco * s.ecco
	omeosq := 1.0 - eccsq
	rteosq := math.Sqrt(omeosq)
	cosio := math.Cos(s.inclo)
	cosio2 := cosio * cosio

	/* ------------------ un-kozai the mean motion ----------------- */
	ak := math.Pow(s.xke/s.no_kozai, X2O3)
	d1 := 0.75 * s.j2 * (3.0*cosio2 - 1.0) / (rteosq * omeosq)
	del := d1 / (ak * ak)
	adel := ak * (1.0 - del*del - del*(1.0/3.0+134.0*del*del/81.0))
	del = d1 / (adel * adel)
	s.no_unkozai = s.no_kozai / (1.0 + del)

	ao := math.Pow(s.xke/s.no_unkozai, X2O3)
	sinio = math.Sin(s.inclo)
	po := ao * omeosq
	con42 := 1.0 - 5.0*cosio2
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

func (s *elsetrec) sgp4init(whichconst string, t *tle, opsmode rune, satn int) {
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

}
