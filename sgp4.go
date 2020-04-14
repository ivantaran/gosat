package gosat

import (
	"errors"
	"math"
)

//Constants
const (
	TwoPi   = math.Pi * 2.0
	Deg2Rad = math.Pi / 180.0
	X2O3    = 2.0 / 3.0
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
	noKozai     float64

	// sgp4fix add new variables from tle
	classification rune
	intldesg       string // TODO: length must be 11
	ephtype        int
	elnum          int // TODO: long originally
	revnum         int // TODO: long originally
	// sgp4fix add unkozai'd variable
	noUnkozai float64
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

type initlVars struct {
	ainv   float64
	ao     float64
	con42  float64
	cosio  float64
	cosio2 float64
	eccsq  float64
	omeosq float64
	posq   float64
	rp     float64
	rteosq float64
	sinio  float64
}

type sgp4ds struct {
	cnodm                   float64
	cosim                   float64
	cosomm                  float64
	day                     float64
	dndt                    float64 // TODO remove this field
	em                      float64
	emsq                    float64
	gam                     float64
	nm                      float64
	mm, nodem, argpm, inclm float64
	rtemsq                  float64
	s1                      float64
	s2                      float64
	s3                      float64
	s4                      float64
	s5                      float64
	s6                      float64
	s7                      float64
	sinim                   float64
	sinomm                  float64
	snodm                   float64
	ss1                     float64
	ss2                     float64
	ss3                     float64
	ss4                     float64
	ss5                     float64
	ss6                     float64
	ss7                     float64
	sz1                     float64
	sz11                    float64
	sz12                    float64
	sz13                    float64
	sz2                     float64
	sz21                    float64
	sz22                    float64
	sz23                    float64
	sz3                     float64
	sz31                    float64
	sz32                    float64
	sz33                    float64
	z1                      float64
	z11                     float64
	z12                     float64
	z13                     float64
	z2                      float64
	z21                     float64
	z22                     float64
	z23                     float64
	z3                      float64
	z31                     float64
	z32                     float64
	z33                     float64
}

type dpperVars struct {
	init                        rune    // read only
	inclo                       float64 // read only
	ep, inclp, nodep, argpp, mp float64
}

//Tle TLE record container
type Tle struct {
	class   byte
	cs1     int
	cs2     int
	design  string
	elsetn  int
	ephtype int
	epoch   float64
	revn    int
	satn    int // catalog number
	title   string
	xargpo  float64
	xbstar  float64
	xecco   float64
	xinclo  float64
	xmo     float64
	xnddot  float64
	xndot   float64
	xno     float64 // mean motion
	xnodeo  float64
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

// inclm, satrec.init, satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo,
// satrec.inclo, 'n', ep, xincp, nodep, argpp, mp,
// double inclo, char init, double &ep, double &inclp, double &nodep, double &argpp, double &mp,
func (s *elsetrec) dpper(vars *dpperVars) {
	// dpper(satrec.e3, satrec.ee2, satrec.peo, satrec.pgho, satrec.pho, satrec.pinco,
	// 	satrec.plo, satrec.se2, satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4,
	// 	satrec.sh2, satrec.sh3, satrec.si2, satrec.si3, satrec.sl2, satrec.sl3,
	// 	satrec.sl4, satrec.t, satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
	// 	satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2, satrec.xl3, satrec.xl4,
	// 	satrec.zmol, satrec.zmos,
	// satrec.operationmode);

	// dpper(satrec.e3, satrec.ee2, satrec.peo, satrec.pgho, satrec.pho, satrec.pinco, satrec.plo,
	// 	satrec.se2, satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4, satrec.sh2, satrec.sh3,
	// 	satrec.si2, satrec.si3, satrec.sl2, satrec.sl3, satrec.sl4, satrec.t, satrec.xgh2,
	// 	satrec.xgh3, satrec.xgh4, satrec.xh2, satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2,
	// 	satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos,
	// satrec.operationmode);

	// static void dpper(double e3, double ee2, double peo, double pgho, double pho, double pinco,
	// 	double plo, double se2, double se3, double sgh2, double sgh3, double sgh4,
	// 	double sh2, double sh3, double si2, double si3, double sl2, double sl3,
	// 	double sl4, double t, double xgh2, double xgh3, double xgh4, double xh2,
	// 	double xh3, double xi2, double xi3, double xl2, double xl3, double xl4,
	// 	double zmol, double zmos,
	// char opsmode) {

	// double alfdp, betdp, cosip, cosop, dalf, dbet, dls, f2, f3, pe, pgh, ph, pinc, pl, sel, ses,
	// sghl, sghs, shll, shs, sil, sinip, sinop, sinzf, sis, sll, sls, xls, xnoh, zf, zm, zel, zes,
	// znl, zns;

	// inclm, satrec.init, satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo,
	// satrec.inclo, 'n', ep, xincp, nodep, argpp, mp,
	// double inclo, char init, double &ep, double &inclp, double &nodep, double &argpp, double &mp,

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
	if vars.init == 'y' { // TODO: else..
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
	if vars.init == 'y' { // TODO: else..
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

	if vars.init == 'n' {
		pe = pe - s.peo
		pinc = pinc - s.pinco
		pl = pl - s.plo
		pgh = pgh - s.pgho
		ph = ph - s.pho
		vars.inclp += pinc
		vars.ep += pe
		sinip := math.Sin(vars.inclp)
		cosip := math.Cos(vars.inclp)

		/* ----------------- apply periodics directly ------------ */
		//  sgp4fix for lyddane choice
		//  strn3 used original inclination - this is technically feasible
		//  gsfc used perturbed inclination - also technically feasible
		//  probably best to readjust the 0.2 limit value and limit discontinuity
		//  0.2 rad = 11.45916 deg
		//  use next line for original strn3 approach and original inclination
		//  if (inclo >= 0.2)
		//  use next line for gsfc version and perturbed inclination
		if vars.inclp >= 0.2 {
			ph = ph / sinip
			pgh = pgh - cosip*ph
			vars.argpp += pgh
			vars.nodep += ph
			vars.mp += pl
		} else {
			/* ---- apply periodics with lyddane modification ---- */
			sinop := math.Sin(vars.nodep)
			cosop := math.Cos(vars.nodep)
			alfdp := sinip * sinop
			betdp := sinip * cosop
			dalf := ph*cosop + pinc*cosip*sinop
			dbet := -ph*sinop + pinc*cosip*cosop
			alfdp = alfdp + dalf
			betdp = betdp + dbet
			vars.nodep = math.Mod(vars.nodep, TwoPi)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (vars.nodep < 0.0) && (s.operationmode == 'a') {
				vars.nodep += TwoPi
			}
			xls := vars.mp + vars.argpp + cosip*vars.nodep
			dls := pl + pgh - pinc*vars.nodep*sinip
			xls = xls + dls
			xnoh := vars.nodep
			vars.nodep = math.Atan2(alfdp, betdp)
			//  sgp4fix for afspc written intrinsic functions
			// nodep used without a trigonometric function ahead
			if (vars.nodep < 0.0) && (s.operationmode == 'a') {
				vars.nodep += TwoPi
			}
			if math.Abs(xnoh-vars.nodep) > math.Pi {
				if vars.nodep < xnoh {
					vars.nodep += TwoPi
				} else {
					vars.nodep -= TwoPi
				}
			}
			vars.mp += pl
			vars.argpp = xls - vars.mp - cosip*vars.nodep
		}
	}
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

func (s *elsetrec) initl(epoch float64) *initlVars {
	var v initlVars
	/* ------------- calculate auxillary epoch quantities ---------- */
	v.eccsq = s.ecco * s.ecco
	v.omeosq = 1.0 - v.eccsq
	v.rteosq = math.Sqrt(v.omeosq)
	v.cosio = math.Cos(s.inclo)
	v.cosio2 = v.cosio * v.cosio

	/* ------------------ un-kozai the mean motion ----------------- */
	ak := math.Pow(s.xke/s.noKozai, X2O3)
	d1 := 0.75 * s.j2 * (3.0*v.cosio2 - 1.0) / (v.rteosq * v.omeosq)
	del := d1 / (ak * ak)
	adel := ak * (1.0 - del*del - del*(1.0/3.0+134.0*del*del/81.0))
	del = d1 / (adel * adel)
	s.noUnkozai = s.noKozai / (1.0 + del)

	v.ao = math.Pow(s.xke/s.noUnkozai, X2O3)
	v.sinio = math.Sin(s.inclo)
	po := v.ao * v.omeosq
	v.con42 = 1.0 - 5.0*v.cosio2
	s.con41 = -v.con42 - v.cosio2 - v.cosio2
	v.ainv = 1.0 / v.ao
	v.posq = po * po
	v.rp = v.ao * (1.0 - s.ecco)
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

	return &v
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
		fallthrough
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
*                           procedure dspace
*
*  this procedure provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
type dspaceVars struct {
	tc,
	argpm,
	dndt,
	em,
	inclm,
	mm,
	nm,
	nodem float64
}

func (vars *dspaceVars) dspace(s *elsetrec) {
	// (
	// 	int irez, double d2201, double d2211, double d3210, double d3222, double d4410,
	// 	double d4422, double d5220, double d5232, double d5421, double d5433,
	// 	double dedt, double del1, double del2, double del3, double didt, double dmdt,
	// 	double dnodt, double domdt, double argpo, double argpdot, double t, double tc,
	// 	double gsto, double xfact, double xlamo, double no, double &atime, double &em,
	// 	double &argpm, double &inclm, double &xli, double &mm, double &xni,
	// 	double &nodem, double &dndt, double &nm) {
	// int iretn, iret;
	// double delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi, g22, g32, g44, g52, g54,
	// fasx2, fasx4, fasx6, rptim, step2, stepn, stepp;

	// s.dspace(s.irez, s.d2201, s.d2211, s.d3210, s.d3222, s.d4410,
	// 	s.d4422, s.d5220, s.d5232, s.d5421, s.d5433, s.dedt,
	// 	s.del1, s.del2, s.del3, s.didt, s.dmdt, s.dnodt,
	// 	s.domdt, s.argpo, s.argpdot, s.t, tc, s.gsto, s.xfact,
	// 	s.xlamo, s.noUnkozai, s.atime, em, argpm, inclm, s.xli, mm,
	// 	s.xni, nodem, dndt, nm)

	const (
		fasx2 = 0.13130908
		fasx4 = 2.8843198
		fasx6 = 0.37448087
		g22   = 5.7686396
		g32   = 0.95240898
		g44   = 1.8014998
		g52   = 1.0508330
		g54   = 4.4108898
		rptim = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
		stepp = 720.0
		stepn = -720.0
		step2 = 259200.0
	)

	/* ----------- calculate deep space resonance effects ----------- */
	vars.dndt = 0.0 // TODO dndt dspaceVars and dssgp4
	theta := math.Mod(s.gsto+vars.tc*rptim, TwoPi)
	vars.em += s.dedt * s.t

	vars.inclm += s.didt * s.t
	vars.argpm += s.domdt * s.t
	vars.nodem += s.dnodt * s.t
	vars.mm += s.dmdt * s.t

	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	//  if (inclm < 0.0)
	// {
	//    inclm = -inclm;
	//    argpm = argpm - pi;
	//    nodem = nodem + pi;
	//  }

	/* - update resonances : numerical (euler-maclaurin) integration - */
	/* ------------------------- epoch restart ----------------------  */
	//   sgp4fix for propagator problems
	//   the following integration works for negative time steps and periods
	//   the specific changes are unknown because the original code was so convoluted

	// sgp4fix take out atime = 0.0 and fix for faster operation
	ft := 0.0
	if s.irez != 0 {
		// sgp4fix streamline check
		if (s.atime == 0.0) || (s.t*s.atime <= 0.0) || (math.Abs(s.t) < math.Abs(s.atime)) {
			s.atime = 0.0
			s.xni = s.noUnkozai // TODO !!!??
			s.xli = s.xlamo
		}
		// sgp4fix move check outside loop
		var delt float64
		if s.t > 0.0 {
			delt = stepp
		} else {
			delt = stepn
		}

		iretn := 381 // added for do loop
		xndt := 0.0
		xldot := 0.0
		xnddt := 0.0
		xomi := 0.0
		x2li := 0.0
		x2omi := 0.0
		for iretn == 381 {
			/* ------------------- dot terms calculated ------------- */
			/* ----------- near - synchronous resonance terms ------- */
			if s.irez != 2 {
				xndt = s.del1*math.Sin(s.xli-fasx2) + s.del2*math.Sin(2.0*(s.xli-fasx4)) +
					s.del3*math.Sin(3.0*(s.xli-fasx6))
				xldot = s.xni + s.xfact
				xnddt = s.del1*math.Cos(s.xli-fasx2) + 2.0*s.del2*math.Cos(2.0*(s.xli-fasx4)) +
					3.0*s.del3*math.Cos(3.0*(s.xli-fasx6))
				xnddt = xnddt * xldot
			} else {
				/* --------- near - half-day resonance terms -------- */
				xomi = s.argpo + s.argpdot*s.atime
				x2omi = xomi + xomi
				x2li = s.xli + s.xli
				xndt = s.d2201*math.Sin(x2omi+s.xli-g22) + s.d2211*math.Sin(s.xli-g22) +
					s.d3210*math.Sin(xomi+s.xli-g32) + s.d3222*math.Sin(-xomi+s.xli-g32) +
					s.d4410*math.Sin(x2omi+x2li-g44) + s.d4422*math.Sin(x2li-g44) +
					s.d5220*math.Sin(xomi+s.xli-g52) + s.d5232*math.Sin(-xomi+s.xli-g52) +
					s.d5421*math.Sin(xomi+x2li-g54) + s.d5433*math.Sin(-xomi+x2li-g54)
				xldot = s.xni + s.xfact
				xnddt = s.d2201*math.Cos(x2omi+s.xli-g22) + s.d2211*math.Cos(s.xli-g22) +
					s.d3210*math.Cos(xomi+s.xli-g32) + s.d3222*math.Cos(-xomi+s.xli-g32) +
					s.d5220*math.Cos(xomi+s.xli-g52) + s.d5232*math.Cos(-xomi+s.xli-g52) +
					2.0*(s.d4410*math.Cos(x2omi+x2li-g44)+s.d4422*math.Cos(x2li-g44)+
						s.d5421*math.Cos(xomi+x2li-g54)+s.d5433*math.Cos(-xomi+x2li-g54))
				xnddt = xnddt * xldot
			}

			/* ----------------------- integrator ------------------- */
			// sgp4fix move end checks to end of routine
			if math.Abs(s.t-s.atime) >= stepp {
				iretn = 381
			} else { // exit here
				ft = s.t - s.atime
				iretn = 0
			}

			if iretn == 381 {
				s.xli += xldot*delt + xndt*step2
				s.xni += xndt*delt + xnddt*step2
				s.atime += delt
			}
		}

		vars.nm = s.xni + xndt*ft + xnddt*ft*ft*0.5
		xl := s.xli + xldot*ft + xndt*ft*ft*0.5
		if s.irez != 1 {
			vars.mm = xl - 2.0*vars.nodem + 2.0*theta
			vars.dndt = vars.nm - s.noUnkozai
		} else {
			vars.mm = xl - vars.nodem - vars.argpm + theta
			vars.dndt = vars.nm - s.noUnkozai
		}
		vars.nm = s.noUnkozai + vars.dndt
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
func (s *elsetrec) dscom(ds *sgp4ds, epoch float64, tc float64) {
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

	ds.nm = s.noUnkozai
	ds.em = s.ecco
	ds.snodm = math.Sin(s.nodeo)
	ds.cnodm = math.Cos(s.nodeo)
	ds.sinomm = math.Sin(s.argpo)
	ds.cosomm = math.Cos(s.argpo)
	ds.sinim = math.Sin(s.inclo)
	ds.cosim = math.Cos(s.inclo)
	ds.emsq = ds.em * ds.em
	betasq := 1.0 - ds.emsq
	ds.rtemsq = math.Sqrt(betasq)

	/* ----------------- initialize lunar solar terms --------------- */
	s.peo = 0.0
	s.pinco = 0.0
	s.plo = 0.0
	s.pgho = 0.0
	s.pho = 0.0
	ds.day = epoch + 18261.5 + tc/1440.0 //TODO remove day from struct, it's local variable
	xnodce := math.Mod(4.5236020-9.2422029e-4*ds.day, TwoPi)
	stem := math.Sin(xnodce)
	ctem := math.Cos(xnodce)
	zcosil := 0.91375164 - 0.03568096*ctem
	zsinil := math.Sqrt(1.0 - zcosil*zcosil)
	zsinhl := 0.089683511 * stem / zsinil
	zcoshl := math.Sqrt(1.0 - zsinhl*zsinhl)
	ds.gam = 5.8351514 + 0.0019443680*ds.day
	zx := 0.39785416 * stem / zsinil
	zy := zcoshl*ctem + 0.91744867*zsinhl*stem
	zx = math.Atan2(zx, zy)
	zx = ds.gam + zx - xnodce
	zcosgl := math.Cos(zx)
	zsingl := math.Sin(zx)

	/* ------------------------- do solar terms --------------------- */
	zcosg := zcosgs
	zsing := zsings
	zcosi := zcosis
	zsini := zsinis
	zcosh := ds.cnodm
	zsinh := ds.snodm
	cc := c1ss
	xnoi := 1.0 / ds.nm

	for lsflg := 1; lsflg <= 2; lsflg++ {
		a1 := zcosg*zcosh + zsing*zcosi*zsinh
		a3 := -zsing*zcosh + zcosg*zcosi*zsinh
		a7 := -zcosg*zsinh + zsing*zcosi*zcosh
		a8 := zsing * zsini
		a9 := zsing*zsinh + zcosg*zcosi*zcosh
		a10 := zcosg * zsini
		a2 := ds.cosim*a7 + ds.sinim*a8
		a4 := ds.cosim*a9 + ds.sinim*a10
		a5 := -ds.sinim*a7 + ds.cosim*a8
		a6 := -ds.sinim*a9 + ds.cosim*a10

		x1 := a1*ds.cosomm + a2*ds.sinomm
		x2 := a3*ds.cosomm + a4*ds.sinomm
		x3 := -a1*ds.sinomm + a2*ds.cosomm
		x4 := -a3*ds.sinomm + a4*ds.cosomm
		x5 := a5 * ds.sinomm
		x6 := a6 * ds.sinomm
		x7 := a5 * ds.cosomm
		x8 := a6 * ds.cosomm

		ds.z31 = 12.0*x1*x1 - 3.0*x3*x3
		ds.z32 = 24.0*x1*x2 - 6.0*x3*x4
		ds.z33 = 12.0*x2*x2 - 3.0*x4*x4
		ds.z1 = 3.0*(a1*a1+a2*a2) + ds.z31*ds.emsq
		ds.z2 = 6.0*(a1*a3+a2*a4) + ds.z32*ds.emsq
		ds.z3 = 3.0*(a3*a3+a4*a4) + ds.z33*ds.emsq
		ds.z11 = -6.0*a1*a5 + ds.emsq*(-24.0*x1*x7-6.0*x3*x5)
		ds.z12 = -6.0*(a1*a6+a3*a5) + ds.emsq*(-24.0*(x2*x7+x1*x8)-6.0*(x3*x6+x4*x5))
		ds.z13 = -6.0*a3*a6 + ds.emsq*(-24.0*x2*x8-6.0*x4*x6)
		ds.z21 = 6.0*a2*a5 + ds.emsq*(24.0*x1*x5-6.0*x3*x7)
		ds.z22 = 6.0*(a4*a5+a2*a6) + ds.emsq*(24.0*(x2*x5+x1*x6)-6.0*(x4*x7+x3*x8))
		ds.z23 = 6.0*a4*a6 + ds.emsq*(24.0*x2*x6-6.0*x4*x8)
		ds.z1 = ds.z1 + ds.z1 + betasq*ds.z31
		ds.z2 = ds.z2 + ds.z2 + betasq*ds.z32
		ds.z3 = ds.z3 + ds.z3 + betasq*ds.z33
		ds.s3 = cc * xnoi
		ds.s2 = -0.5 * ds.s3 / ds.rtemsq
		ds.s4 = ds.s3 * ds.rtemsq
		ds.s1 = -15.0 * ds.em * ds.s4
		ds.s5 = x1*x3 + x2*x4
		ds.s6 = x2*x3 + x1*x4
		ds.s7 = x2*x4 - x1*x3

		/* ----------------------- do lunar terms ------------------- */
		if lsflg == 1 {
			ds.ss1 = ds.s1
			ds.ss2 = ds.s2
			ds.ss3 = ds.s3
			ds.ss4 = ds.s4
			ds.ss5 = ds.s5
			ds.ss6 = ds.s6
			ds.ss7 = ds.s7
			ds.sz1 = ds.z1
			ds.sz2 = ds.z2
			ds.sz3 = ds.z3
			ds.sz11 = ds.z11
			ds.sz12 = ds.z12
			ds.sz13 = ds.z13
			ds.sz21 = ds.z21
			ds.sz22 = ds.z22
			ds.sz23 = ds.z23
			ds.sz31 = ds.z31
			ds.sz32 = ds.z32
			ds.sz33 = ds.z33
			zcosg = zcosgl
			zsing = zsingl
			zcosi = zcosil
			zsini = zsinil
			zcosh = zcoshl*ds.cnodm + zsinhl*ds.snodm
			zsinh = ds.snodm*zcoshl - ds.cnodm*zsinhl
			cc = c1l
		}
	}

	s.zmol = math.Mod(4.7199672+0.22997150*ds.day-ds.gam, TwoPi)
	s.zmos = math.Mod(6.2565837+0.017201977*ds.day, TwoPi)

	/* ------------------------ do solar terms ---------------------- */
	s.se2 = 2.0 * ds.ss1 * ds.ss6
	s.se3 = 2.0 * ds.ss1 * ds.ss7
	s.si2 = 2.0 * ds.ss2 * ds.sz12
	s.si3 = 2.0 * ds.ss2 * (ds.sz13 - ds.sz11)
	s.sl2 = -2.0 * ds.ss3 * ds.sz2
	s.sl3 = -2.0 * ds.ss3 * (ds.sz3 - ds.sz1)
	s.sl4 = -2.0 * ds.ss3 * (-21.0 - 9.0*ds.emsq) * zes
	s.sgh2 = 2.0 * ds.ss4 * ds.sz32
	s.sgh3 = 2.0 * ds.ss4 * (ds.sz33 - ds.sz31)
	s.sgh4 = -18.0 * ds.ss4 * zes
	s.sh2 = -2.0 * ds.ss2 * ds.sz22
	s.sh3 = -2.0 * ds.ss2 * (ds.sz23 - ds.sz21)

	/* ------------------------ do lunar terms ---------------------- */
	s.ee2 = 2.0 * ds.s1 * ds.s6
	s.e3 = 2.0 * ds.s1 * ds.s7
	s.xi2 = 2.0 * ds.s2 * ds.z12
	s.xi3 = 2.0 * ds.s2 * (ds.z13 - ds.z11)
	s.xl2 = -2.0 * ds.s3 * ds.z2
	s.xl3 = -2.0 * ds.s3 * (ds.z3 - ds.z1)
	s.xl4 = -2.0 * ds.s3 * (-21.0 - 9.0*ds.emsq) * zel
	s.xgh2 = 2.0 * ds.s4 * ds.z32
	s.xgh3 = 2.0 * ds.s4 * (ds.z33 - ds.z31)
	s.xgh4 = -18.0 * ds.s4 * zel
	s.xh2 = -2.0 * ds.s2 * ds.z22
	s.xh3 = -2.0 * ds.s2 * (ds.z23 - ds.z21)
}

/*-----------------------------------------------------------------------------
*
*                           procedure dsinit
*
*  this procedure provides deep space contributions to mean motion dot due
*    to geopotential resonance with half day and one day orbits.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    xke         - reciprocal of tumin
*    cosim, sinim-
*    emsq        - eccentricity squared
*    argpo       - argument of perigee
*    s1, s2, s3, s4, s5      -
*    ss1, ss2, ss3, ss4, ss5 -
*    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
*    t           - time
*    tc          -
*    gsto        - greenwich sidereal time                   rad
*    mo          - mean anomaly
*    mdot        - mean anomaly dot (rate)
*    no          - mean motion
*    nodeo       - right ascension of ascending node
*    nodedot     - right ascension of ascending node dot (rate)
*    xpidot      -
*    z1, z3, z11, z13, z21, z23, z31, z33 -
*    eccm        - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    xn          - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right ascension of ascending node
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    atime       -
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
*    dedt        -
*    didt        -
*    dmdt        -
*    dndt        -
*    dnodt       -
*    domdt       -
*    del1, del2, del3        -
*    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
*    theta       -
*    xfact       -
*    xlamo       -
*    xli         -
*    xni
*
*  locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
*    sini2       -
*    temp        -
*    temp1       -
*    theta       -
*    xno2        -
*
*  coupling      :
*    getgravconst- no longer used
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/
// type dsinitVars struct {
// 	d2201 float64
// 	d2211 float64
// 	d3210 float64
// 	d3222 float64
// 	d4410 float64
// 	d4422 float64
// 	d5220 float64
// 	d5232 float64
// 	d5421 float64
// 	d5433 float64
// 	dedt  float64
// 	del1  float64
// 	del2  float64
// 	del3  float64
// 	didt  float64
// 	dmdt  float64
// 	dndt  float64
// 	dnodt float64
// 	domdt float64
// 	xfact float64
// 	xlamo float64
// 	xli   float64
// 	xni   float64
// }

func (s *elsetrec) dsinit(ds *sgp4ds, iv *initlVars, tc, xpidot float64) {
	// dsinit(satrec.xke, cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3,
	// 	ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec.t, tc,
	// 	satrec.gsto, satrec.mo, satrec.mdot, satrec.no_unkozai, satrec.nodeo,
	// 	satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33, satrec.ecco, eccsq,
	// 	em, argpm, inclm, mm, nm, nodem, satrec.irez, satrec.atime, satrec.d2201,
	// 	satrec.d2211, satrec.d3210, satrec.d3222, satrec.d4410, satrec.d4422,
	// 	satrec.d5220, satrec.d5232, satrec.d5421, satrec.d5433, satrec.dedt, satrec.didt,
	// 	satrec.dmdt, dndt, satrec.dnodt, satrec.domdt, satrec.del1, satrec.del2,
	// 	satrec.del3, satrec.xfact, satrec.xlamo, satrec.xli, satrec.xni);

	// // sgp4fix just send in xke as a constant and eliminate getgravconst call
	// // gravconsttype whichconst,
	// double xke, double cosim, double emsq, double argpo, double s1, double s2, double s3, double s4,
	// double s5, double sinim, double ss1, double ss2, double ss3, double ss4, double ss5, double sz1,
	// double sz3, double sz11, double sz13, double sz21, double sz23, double sz31, double sz33,
	// double t, double tc, double gsto, double mo, double mdot, double no, double nodeo,
	// double nodedot, double xpidot, double z1, double z3, double z11, double z13, double z21,
	// double z23, double z31, double z33, double ecco, double eccsq, double &em, double &argpm,
	// double &inclm, double &mm, double &nm, double &nodem, int &irez, double &atime, double &d2201,
	// double &d2211, double &d3210, double &d3222, double &d4410, double &d4422, double &d5220,
	// double &d5232, double &d5421, double &d5433, double &dedt, double &didt, double &dmdt,
	// double &dndt, double &dnodt, double &domdt, double &del1, double &del2, double &del3,
	// double &xfact, double &xlamo, double &xli, double &xni) {
	/* --------------------- local variables ------------------------ */

	// double ainv2, aonv = 0.0, cosisq, eoc, f220, f221, f311, f321, f322, f330, f441, f442, f522,
	//               f523, f542, f543, g200, g201, g211, g300, g310, g322, g410, g422, g520, g521,
	//               g532, g533, ses, sgs, sghl, sghs, shs, shll, sis, sini2, sls, temp, temp1, theta,
	//               xno2, q22, q31, q33, root22, root44, root54, rptim, root32, root52, x2o3, znl,
	//               emo, zns, emsqo;

	aonv := 0.0
	const (
		q22    = 1.7891679e-6
		q31    = 2.1460748e-6
		q33    = 2.2123015e-7
		root22 = 1.7891679e-6
		root44 = 7.3636953e-9
		root54 = 2.1765803e-9
		rptim  = 4.37526908801129966e-3 // this equates to 7.29211514668855e-5 rad/sec
		root32 = 3.7393792e-7
		root52 = 1.1428639e-7
		znl    = 1.5835218e-4
		zns    = 1.19459e-5
	)
	// sgp4fix identify constants and allow alternate values
	// just xke is used here so pass it in rather than have multiple calls
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

	/* -------------------- deep space initialization ------------ */
	s.irez = 0 // TODO: unwrap if's
	if (ds.nm < 0.0052359877) && (ds.nm > 0.0034906585) {
		s.irez = 1
	}
	if (ds.nm >= 8.26e-3) && (ds.nm <= 9.24e-3) && (ds.em >= 0.5) {
		s.irez = 2
	}

	/* ------------------------ do solar terms ------------------- */
	ses := ds.ss1 * zns * ds.ss5
	sis := ds.ss2 * zns * (ds.sz11 + ds.sz13)
	sls := -zns * ds.ss3 * (ds.sz1 + ds.sz3 - 14.0 - 6.0*ds.emsq)
	sghs := ds.ss4 * zns * (ds.sz31 + ds.sz33 - 6.0)
	shs := -zns * ds.ss2 * (ds.sz21 + ds.sz23)
	// sgp4fix for 180 deg incl
	if (ds.inclm < 5.2359877e-2) || (ds.inclm > math.Pi-5.2359877e-2) {
		shs = 0.0
	}
	if ds.sinim != 0.0 {
		shs = shs / ds.sinim
	}
	sgs := sghs - ds.cosim*shs

	/* ------------------------- do lunar terms ------------------ */
	s.dedt = ses + ds.s1*znl*ds.s5
	s.didt = sis + ds.s2*znl*(ds.z11+ds.z13)
	s.dmdt = sls - znl*ds.s3*(ds.z1+ds.z3-14.0-6.0*ds.emsq)
	sghl := ds.s4 * znl * (ds.z31 + ds.z33 - 6.0)
	shll := -znl * ds.s2 * (ds.z21 + ds.z23)
	// sgp4fix for 180 deg incl
	if (ds.inclm < 5.2359877e-2) || (ds.inclm > math.Pi-5.2359877e-2) {
		shll = 0.0
	}
	s.domdt = sgs + sghl
	s.dnodt = shs
	if ds.sinim != 0.0 { // TODO: comparing float with zero
		s.domdt -= ds.cosim / ds.sinim * shll
		s.dnodt += shll / ds.sinim
	}

	/* ----------- calculate deep space resonance effects -------- */
	ds.dndt = 0.0
	theta := math.Mod(s.gsto+tc*rptim, TwoPi)
	ds.em += s.dedt * s.t
	ds.inclm += s.didt * s.t
	ds.argpm += s.domdt * s.t
	ds.nodem += s.dnodt * s.t
	ds.mm += s.dmdt * s.t
	//   sgp4fix for negative inclinations
	//   the following if statement should be commented out
	// if (inclm < 0.0)
	//  {
	//    inclm  = -inclm;
	//    argpm  = argpm - pi;
	//    nodem = nodem + pi;
	//  }

	/* -------------- initialize the resonance terms ------------- */
	if s.irez != 0 {
		aonv = math.Pow(ds.nm/s.xke, X2O3)

		/* ---------- geopotential resonance for 12 hour orbits ------ */
		if s.irez == 2 {
			cosisq := ds.cosim * ds.cosim
			emo := ds.em
			ds.em = s.ecco
			emsqo := ds.emsq
			ds.emsq = iv.eccsq
			eoc := ds.em * ds.emsq
			g201 := -0.306 - (ds.em-0.64)*0.440

			var g211, g310, g322, g410, g422, g520, g521, g532, g533 float64

			if ds.em <= 0.65 {
				g211 = 3.616 - 13.2470*ds.em + 16.2900*ds.emsq
				g310 = -19.302 + 117.3900*ds.em - 228.4190*ds.emsq + 156.5910*eoc
				g322 = -18.9068 + 109.7927*ds.em - 214.6334*ds.emsq + 146.5816*eoc
				g410 = -41.122 + 242.6940*ds.em - 471.0940*ds.emsq + 313.9530*eoc
				g422 = -146.407 + 841.8800*ds.em - 1629.014*ds.emsq + 1083.4350*eoc
				g520 = -532.114 + 3017.977*ds.em - 5740.032*ds.emsq + 3708.2760*eoc
			} else {
				g211 = -72.099 + 331.819*ds.em - 508.738*ds.emsq + 266.724*eoc
				g310 = -346.844 + 1582.851*ds.em - 2415.925*ds.emsq + 1246.113*eoc
				g322 = -342.585 + 1554.908*ds.em - 2366.899*ds.emsq + 1215.972*eoc
				g410 = -1052.797 + 4758.686*ds.em - 7193.992*ds.emsq + 3651.957*eoc
				g422 = -3581.690 + 16178.110*ds.em - 24462.770*ds.emsq + 12422.520*eoc
				if ds.em > 0.715 {
					g520 = -5149.66 + 29936.92*ds.em - 54087.36*ds.emsq + 31324.56*eoc
				} else {
					g520 = 1464.74 - 4664.75*ds.em + 3763.64*ds.emsq
				}
			}
			if ds.em < 0.7 {
				g533 = -919.22770 + 4988.6100*ds.em - 9064.7700*ds.emsq + 5542.21*eoc
				g521 = -822.71072 + 4568.6173*ds.em - 8491.4146*ds.emsq + 5337.524*eoc
				g532 = -853.66600 + 4690.2500*ds.em - 8624.7700*ds.emsq + 5341.4*eoc
			} else {
				g533 = -37995.780 + 161616.52*ds.em - 229838.20*ds.emsq + 109377.94*eoc
				g521 = -51752.104 + 218913.95*ds.em - 309468.16*ds.emsq + 146349.42*eoc
				g532 = -40023.880 + 170470.89*ds.em - 242699.48*ds.emsq + 115605.82*eoc
			}

			sini2 := ds.sinim * ds.sinim
			f220 := 0.75 * (1.0 + 2.0*ds.cosim + cosisq)
			f221 := 1.5 * sini2
			f321 := 1.875 * ds.sinim * (1.0 - 2.0*ds.cosim - 3.0*cosisq)
			f322 := -1.875 * ds.sinim * (1.0 + 2.0*ds.cosim - 3.0*cosisq)
			f441 := 35.0 * sini2 * f220
			f442 := 39.3750 * sini2 * sini2
			f522 := 9.84375 * ds.sinim * (sini2*(1.0-2.0*ds.cosim-5.0*cosisq) +
				0.33333333*(-2.0+4.0*ds.cosim+6.0*cosisq))
			f523 := ds.sinim * (4.92187512*sini2*(-2.0-4.0*ds.cosim+10.0*cosisq) +
				6.56250012*(1.0+2.0*ds.cosim-3.0*cosisq))
			f542 := 29.53125 * ds.sinim *
				(2.0 - 8.0*ds.cosim + cosisq*(-12.0+8.0*ds.cosim+10.0*cosisq))
			f543 := 29.53125 * ds.sinim *
				(-2.0 - 8.0*ds.cosim + cosisq*(12.0+8.0*ds.cosim-10.0*cosisq))
			xno2 := ds.nm * ds.nm
			ainv2 := aonv * aonv
			temp1 := 3.0 * xno2 * ainv2
			temp := temp1 * root22
			s.d2201 = temp * f220 * g201
			s.d2211 = temp * f221 * g211
			temp1 = temp1 * aonv
			temp = temp1 * root32
			s.d3210 = temp * f321 * g310
			s.d3222 = temp * f322 * g322
			temp1 = temp1 * aonv
			temp = 2.0 * temp1 * root44
			s.d4410 = temp * f441 * g410
			s.d4422 = temp * f442 * g422
			temp1 = temp1 * aonv
			temp = temp1 * root52
			s.d5220 = temp * f522 * g520
			s.d5232 = temp * f523 * g532
			temp = 2.0 * temp1 * root54
			s.d5421 = temp * f542 * g521
			s.d5433 = temp * f543 * g533
			s.xlamo = math.Mod(s.mo+s.nodeo+s.nodeo-theta-theta, TwoPi)
			s.xfact = s.mdot + s.dmdt + 2.0*(s.nodedot+s.dnodt-rptim) - s.noUnkozai
			ds.em = emo
			ds.emsq = emsqo
		}

		/* ---------------- synchronous resonance terms -------------- */
		if s.irez == 1 {
			g200 := 1.0 + ds.emsq*(-2.5+0.8125*ds.emsq)
			g310 := 1.0 + 2.0*ds.emsq
			g300 := 1.0 + ds.emsq*(-6.0+6.60937*ds.emsq)
			f220 := 0.75 * (1.0 + ds.cosim) * (1.0 + ds.cosim)
			f311 := 0.9375*ds.sinim*ds.sinim*(1.0+3.0*ds.cosim) - 0.75*(1.0+ds.cosim)
			f330 := 1.0 + ds.cosim
			f330 = 1.875 * f330 * f330 * f330
			s.del1 = 3.0 * ds.nm * ds.nm * aonv * aonv
			s.del2 = 2.0 * s.del1 * f220 * g200 * q22
			s.del3 = 3.0 * s.del1 * f330 * g300 * q33 * aonv
			s.del1 = s.del1 * f311 * g310 * q31 * aonv
			s.xlamo = math.Mod(s.mo+s.nodeo+s.argpo-theta, TwoPi)
			s.xfact = s.mdot + xpidot - rptim + s.dmdt + s.domdt + s.dnodt - s.noUnkozai
		}

		/* ------------ for sgp4, initialize the integrator ---------- */
		s.xli = s.xlamo
		s.xni = s.noUnkozai
		s.atime = 0.0
		ds.nm = s.noUnkozai + ds.dndt
	}
}

func (s *elsetrec) sgp4init(whichconst string, t *Tle, opsmode rune) error {
	const temp4 = 1.5e-12

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
	s.satnum = t.satn

	s.bstar = t.xbstar
	s.ndot = t.xndot
	s.nddot = t.xnddot
	s.ecco = t.xecco
	s.argpo = t.xargpo
	s.inclo = t.xinclo
	s.mo = t.xmo
	s.noKozai = t.xno
	s.nodeo = t.xnodeo

	s.am = 0.0
	s.em = 0.0
	s.im = 0.0
	s.Om = 0.0
	s.om = 0.0
	s.mm = 0.0
	s.nm = 0.0

	ss := 78.0/s.radiusearthkm + 1.0
	qzms2ttemp := (120.0 - 78.0) / s.radiusearthkm
	qzms2t := qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp

	s.init = 'y'
	s.t = 0.0

	iv := s.initl(t.epoch)

	s.a = math.Pow(s.noUnkozai*s.tumin, -X2O3)
	s.alta = s.a*(1.0+s.ecco) - 1.0
	s.altp = s.a*(1.0-s.ecco) - 1.0
	s.error = 0

	if (iv.omeosq >= 0.0) || (s.noUnkozai >= 0.0) {
		s.isimp = 0
		if iv.rp < (220.0/s.radiusearthkm + 1.0) {
			s.isimp = 1
		}
		sfour := ss
		qzms24 := qzms2t
		perige := (iv.rp - 1.0) * s.radiusearthkm

		/* - for perigees below 156 km, s and qoms2t are altered - */
		if perige < 156.0 {
			sfour = perige - 78.0
			if perige < 98.0 {
				sfour = 20.0
			}
			qzms24temp := (120.0 - sfour) / s.radiusearthkm
			qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp
			sfour = sfour/s.radiusearthkm + 1.0
		}
		pinvsq := 1.0 / iv.posq

		tsi := 1.0 / (iv.ao - sfour)
		s.eta = iv.ao * s.ecco * tsi
		etasq := s.eta * s.eta
		eeta := s.ecco * s.eta
		psisq := math.Abs(1.0 - etasq)
		coef := qzms24 * math.Pow(tsi, 4.0)
		coef1 := coef / math.Pow(psisq, 3.5)
		cc2 := coef1 * s.noUnkozai * (iv.ao*(1.0+1.5*etasq+eeta*(4.0+etasq)) +
			0.375*s.j2*tsi/psisq*s.con41*(8.0+3.0*etasq*(8.0+etasq)))
		s.cc1 = s.bstar * cc2
		cc3 := 0.0
		if s.ecco > 1.0e-4 {
			cc3 = -2.0 * coef * tsi * s.j3oj2 * s.noUnkozai * iv.sinio / s.ecco
		}
		s.x1mth2 = 1.0 - iv.cosio2
		s.cc4 = 2.0 * s.noUnkozai * coef1 * iv.ao * iv.omeosq * (s.eta*(2.0+0.5*etasq) +
			s.ecco*(0.5+2.0*etasq) -
			s.j2*tsi/(iv.ao*psisq)*(-3.0*s.con41*(1.0-2.0*eeta+etasq*(1.5-0.5*eeta))+
				0.75*s.x1mth2*(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*s.argpo)))
		s.cc5 = 2.0 * coef1 * iv.ao * iv.omeosq * (1.0 + 2.75*(etasq+eeta) + eeta*etasq)
		cosio4 := iv.cosio2 * iv.cosio2
		temp1 := 1.5 * s.j2 * pinvsq * s.noUnkozai
		temp2 := 0.5 * temp1 * s.j2 * pinvsq
		temp3 := -0.46875 * s.j4 * pinvsq * pinvsq * s.noUnkozai
		s.mdot = s.noUnkozai + 0.5*temp1*iv.rteosq*s.con41 +
			0.0625*temp2*iv.rteosq*(13.0-78.0*iv.cosio2+137.0*cosio4)
		s.argpdot = -0.5*temp1*iv.con42 + 0.0625*temp2*(7.0-114.0*iv.cosio2+395.0*cosio4) +
			temp3*(3.0-36.0*iv.cosio2+49.0*cosio4)
		xhdot1 := -temp1 * iv.cosio
		s.nodedot = xhdot1 + (0.5*temp2*(4.0-19.0*iv.cosio2)+2.0*temp3*(3.0-7.0*iv.cosio2))*iv.cosio
		xpidot := s.argpdot + s.nodedot
		s.omgcof = s.bstar * cc3 * math.Cos(s.argpo)
		s.xmcof = 0.0
		if s.ecco > 1.0e-4 {
			s.xmcof = -X2O3 * coef * s.bstar / eeta
		}
		s.nodecf = 3.5 * iv.omeosq * xhdot1 * s.cc1
		s.t2cof = 1.5 * s.cc1

		if math.Abs(iv.cosio+1.0) > 1.5e-12 {
			s.xlcof = -0.25 * s.j3oj2 * iv.sinio * (3.0 + 5.0*iv.cosio) / (1.0 + iv.cosio)
		} else {
			s.xlcof = -0.25 * s.j3oj2 * iv.sinio * (3.0 + 5.0*iv.cosio) / temp4
		}
		s.aycof = -0.5 * s.j3oj2 * iv.sinio
		delmotemp := 1.0 + s.eta*math.Cos(s.mo)
		s.delmo = delmotemp * delmotemp * delmotemp
		s.sinmao = math.Sin(s.mo)
		s.x7thm1 = 7.0*iv.cosio2 - 1.0

		/* --------------- deep space initialization ------------- */
		if (TwoPi / s.noUnkozai) >= 225.0 {
			s.method = 'd'
			s.isimp = 1
			tc := 0.0
			var ds sgp4ds
			ds.inclm = s.inclo
			s.dscom(&ds, t.epoch, tc)

			var dpperVars dpperVars
			dpperVars.init = s.init
			dpperVars.inclo = ds.inclm
			dpperVars.ep = s.ecco
			dpperVars.inclp = s.inclo
			dpperVars.nodep = s.nodeo
			dpperVars.argpp = s.argpo
			dpperVars.mp = s.mo

			// dpper(satrec.e3, satrec.ee2, satrec.peo, satrec.pgho, satrec.pho, satrec.pinco,
			// 	satrec.plo, satrec.se2, satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4,
			// 	satrec.sh2, satrec.sh3, satrec.si2, satrec.si3, satrec.sl2, satrec.sl3,
			// 	satrec.sl4, satrec.t, satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2,
			// 	satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2, satrec.xl3, satrec.xl4,
			// 	satrec.zmol, satrec.zmos, inclm, satrec.init, satrec.ecco, satrec.inclo,
			// 	satrec.nodeo, satrec.argpo, satrec.mo, satrec.operationmode);

			// static void dpper(double e3, double ee2, double peo, double pgho, double pho, double pinco,
			// 	double plo, double se2, double se3, double sgh2, double sgh3, double sgh4,
			// 	double sh2, double sh3, double si2, double si3, double sl2, double sl3,
			// 	double sl4, double t, double xgh2, double xgh3, double xgh4, double xh2,
			// 	double xh3, double xi2, double xi3, double xl2, double xl3, double xl4,
			// 	double zmol, double zmos, double inclo, char init, double &ep, double &inclp,
			// 	double &nodep, double &argpp, double &mp, char opsmode) {

			s.dpper(&dpperVars)
			s.ecco = dpperVars.ep
			s.inclo = dpperVars.inclp
			s.nodeo = dpperVars.nodep
			s.argpo = dpperVars.argpp
			s.mo = dpperVars.mp

			ds.argpm = 0.0
			ds.nodem = 0.0
			ds.mm = 0.0
			s.dsinit(&ds, iv, tc, xpidot)
		}

		/* ----------- set variables if not deep space ----------- */
		if s.isimp != 1 {
			cc1sq := s.cc1 * s.cc1
			s.d2 = 4.0 * iv.ao * tsi * cc1sq
			temp := s.d2 * tsi * s.cc1 / 3.0
			s.d3 = (17.0*iv.ao + sfour) * temp
			s.d4 = 0.5 * temp * iv.ao * tsi * (221.0*iv.ao + 31.0*sfour) * s.cc1
			s.t3cof = s.d2 + 2.0*cc1sq
			s.t4cof = 0.25 * (3.0*s.d3 + s.cc1*(12.0*s.d2+10.0*cc1sq))
			s.t5cof = 0.2 * (3.0*s.d4 + 12.0*s.cc1*s.d3 + 6.0*s.d2*s.d2 + 15.0*cc1sq*(2.0*s.d2+cc1sq))
		}
	}

	/* finally propogate to zero epoch to initialize all others. */
	// sgp4fix take out check to let satellites process until they are actually below earth surface
	//       if(satrec.error == 0)
	var r [6]float64
	err := s.sgp4(0.0, r[0:3], r[3:])

	s.init = 'n'

	err = s.sgp4(2880.0, r[0:3], r[3:])

	return err
}

/*-----------------------------------------------------------------------------
*
*                             procedure sgp4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the
*    methodology from the aiaa paper (2006) describing the history and
*    development of the code.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satrec	 - initialised structure from sgp4init() call.
*    tsince	 - time since epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    am          -
*    axnl, aynl        -
*    betal       -
*    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
*    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
*    cosisq  , cossu   , sinsu   , cosu    , sinu
*    delm        -
*    delomg      -
*    dndt        -
*    eccm        -
*    emsq        -
*    ecose       -
*    el2         -
*    eo1         -
*    eccp        -
*    esine       -
*    argpm       -
*    argpp       -
*    omgadf      -c
*    pl          -
*    r           -
*    rtemsq      -
*    rdotl       -
*    rl          -
*    rvdot       -
*    rvdotl      -
*    su          -
*    t2  , t3   , t4    , tc
*    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
*    u   , ux   , uy    , uz     , vx     , vy     , vz
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - right asc of ascending node
*    xinc        -
*    xincp       -
*    xl          -
*    xlm         -
*    mp          -
*    xmdf        -
*    xmx         -
*    xmy         -
*    nodedf      -
*    xnode       -
*    nodep       -
*    np          -
*
*  coupling      :
*    getgravconst- no longer used. Variables are conatined within satrec
*    dpper
*    dpspace
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
----------------------------------------------------------------------------*/

func (s *elsetrec) sgp4(tsince float64, r []float64, v []float64) error {

	// double am, axnl, aynl, betal, cosim, cnod, cos2u, coseo1, cosi, cosip, cosisq, cossu, cosu,
	//     delm, delomg, em, emsq, ecose, el2, eo1, ep, esine, argpm, argpp, argpdf, pl,
	//     mrt = 0.0, mvt, rdotl, rl, rvdot, rvdotl, sinim, sin2u, sineo1, sini, sinip, sinsu, sinu,
	//     snod, su, t2, t3, t4, tem5, temp, temp1, temp2, tempa, tempe, templ, u, ux, uy, uz, vx, vy,
	//     vz, inclm, mm, nm, nodem, xinc, xincp, xl, xlm, mp, xmdf, xmx, xmy, nodedf, xnode, nodep,
	//     tc, dndt, twopi, x2o3, vkmpersec, delmtemp;
	// int ktr;

	/* ------------------ set mathematical constants --------------- */
	// sgp4fix divisor for divide by zero check on inclination
	// the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
	// 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
	const temp4 = 1.5e-12

	var dspaceVars dspaceVars

	// sgp4fix identify constants and allow alternate values
	// getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
	vkmpersec := s.radiusearthkm * s.xke / 60.0

	/* --------------------- clear sgp4 error flag ----------------- */
	s.t = tsince
	s.error = 0

	/* ------- update for secular gravity and atmospheric drag ----- */
	xmdf := s.mo + s.mdot*s.t
	argpdf := s.argpo + s.argpdot*s.t
	nodedf := s.nodeo + s.nodedot*s.t
	dspaceVars.argpm = argpdf
	dspaceVars.mm = xmdf
	t2 := s.t * s.t
	dspaceVars.nodem = nodedf + s.nodecf*t2
	tempa := 1.0 - s.cc1*s.t
	tempe := s.bstar * s.cc4 * s.t
	templ := s.t2cof * t2

	if s.isimp != 1 {
		delomg := s.omgcof * s.t
		// sgp4fix use mutliply for speed instead of pow
		delmtemp := 1.0 + s.eta*math.Cos(xmdf)
		delm := s.xmcof * (delmtemp*delmtemp*delmtemp - s.delmo)
		temp := delomg + delm
		dspaceVars.mm = xmdf + temp
		dspaceVars.argpm = argpdf - temp
		t3 := t2 * s.t
		t4 := t3 * s.t
		tempa = tempa - s.d2*t2 - s.d3*t3 - s.d4*t4
		tempe = tempe + s.bstar*s.cc5*(math.Sin(dspaceVars.mm)-s.sinmao)
		templ = templ + s.t3cof*t3 + t4*(s.t4cof+s.t*s.t5cof)
	}

	dspaceVars.nm = s.noUnkozai
	dspaceVars.em = s.ecco
	dspaceVars.inclm = s.inclo
	if s.method == 'd' {
		dspaceVars.tc = s.t
		dspaceVars.dspace(s)
	}

	if dspaceVars.nm <= 0.0 {
		s.error = 2
		// sgp4fix add return
		return errors.New("nm <= 0.0")
	}
	am := math.Pow(s.xke/dspaceVars.nm, X2O3) * tempa * tempa
	dspaceVars.nm = s.xke / math.Pow(am, 1.5)
	dspaceVars.em -= tempe

	// fix tolerance for error recognition
	// sgp4fix am is fixed from the previous nm check
	if (dspaceVars.em >= 1.0) || (dspaceVars.em < -0.001) /* || (am < 0.95)*/ {
		s.error = 1
		// sgp4fix to return if there is an error in eccentricity
		return errors.New("(em >= 1.0) || (em < -0.001)")
	}
	// sgp4fix fix tolerance to avoid a divide by zero
	if dspaceVars.em < 1.0e-6 {
		dspaceVars.em = 1.0e-6
	}
	dspaceVars.mm += s.noUnkozai * templ
	xlm := dspaceVars.mm + dspaceVars.argpm + dspaceVars.nodem
	emsq := dspaceVars.em * dspaceVars.em
	temp := 1.0 - emsq

	dspaceVars.nodem = math.Mod(dspaceVars.nodem, TwoPi)
	dspaceVars.argpm = math.Mod(dspaceVars.argpm, TwoPi)
	xlm = math.Mod(xlm, TwoPi)
	dspaceVars.mm = math.Mod(xlm-dspaceVars.argpm-dspaceVars.nodem, TwoPi)

	// sgp4fix recover singly averaged mean elements
	s.am = am
	s.em = dspaceVars.em
	s.im = dspaceVars.inclm
	s.Om = dspaceVars.nodem
	s.om = dspaceVars.argpm
	s.mm = dspaceVars.mm
	s.nm = dspaceVars.nm

	/* ----------------- compute extra mean quantities ------------- */
	sinim := math.Sin(dspaceVars.inclm)
	cosim := math.Cos(dspaceVars.inclm)

	/* -------------------- add lunar-solar periodics -------------- */
	ep := dspaceVars.em
	xincp := dspaceVars.inclm
	argpp := dspaceVars.argpm
	nodep := dspaceVars.nodem
	mp := dspaceVars.mm
	sinip := sinim
	cosip := cosim

	var dpperVars dpperVars
	dpperVars.init = 'n'
	dpperVars.inclo = s.inclo
	dpperVars.ep = ep
	dpperVars.inclp = xincp
	dpperVars.nodep = nodep
	dpperVars.argpp = argpp
	dpperVars.mp = mp

	if s.method == 'd' {
		// inclm, satrec.init, satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo,
		// satrec.inclo, 'n', ep, xincp, nodep, argpp, mp,

		// double inclo, char init, double &ep, double &inclp, double &nodep, double &argpp, double &mp,
		// func (s *elsetrec) dpper(vars *dpperVars) {

		s.dpper(&dpperVars)
		ep = dpperVars.ep
		xincp = dpperVars.inclp
		nodep = dpperVars.nodep
		argpp = dpperVars.argpp
		mp = dpperVars.mp

		// dpper(s.e3, s.ee2, s.peo, s.pgho, s.pho, s.pinco, s.plo,
		// 	s.se2, s.se3, s.sgh2, s.sgh3, s.sgh4, s.sh2, s.sh3,
		// 	s.si2, s.si3, s.sl2, s.sl3, s.sl4, s.t, s.xgh2,
		// 	s.xgh3, s.xgh4, s.xh2, s.xh3, s.xi2, s.xi3, s.xl2,
		// 	s.xl3, s.xl4, s.zmol, s.zmos, s.inclo, 'n', ep, xincp, nodep,
		// 	argpp, mp, s.operationmode)
		if xincp < 0.0 {
			xincp = -xincp
			nodep = nodep + math.Pi
			argpp = argpp - math.Pi
		}
		if (ep < 0.0) || (ep > 1.0) {
			s.error = 3
			// sgp4fix add return
			return errors.New("(ep < 0.0) || (ep > 1.0)")
		}
	}

	/* -------------------- long period periodics ------------------ */
	if s.method == 'd' {
		sinip = math.Sin(xincp)
		cosip = math.Cos(xincp)
		s.aycof = -0.5 * s.j3oj2 * sinip
		// sgp4fix for divide by zero for xincp = 180 deg
		if math.Abs(cosip+1.0) > 1.5e-12 {
			s.xlcof = -0.25 * s.j3oj2 * sinip * (3.0 + 5.0*cosip) / (1.0 + cosip)
		} else {
			s.xlcof = -0.25 * s.j3oj2 * sinip * (3.0 + 5.0*cosip) / temp4
		}
	}
	axnl := ep * math.Cos(argpp)
	temp = 1.0 / (am * (1.0 - ep*ep))
	aynl := ep*math.Sin(argpp) + temp*s.aycof
	xl := mp + argpp + nodep + temp*s.xlcof*axnl

	/* --------------------- solve kepler's equation --------------- */
	u := math.Mod(xl-nodep, TwoPi)
	eo1 := u
	tem5 := 9999.9
	ktr := 1
	//   sgp4fix for kepler iteration
	//   the following iteration needs better limits on corrections
	sineo1 := 0.0
	coseo1 := 0.0
	for (math.Abs(tem5) >= 1.0e-12) && (ktr <= 10) {
		sineo1 = math.Sin(eo1)
		coseo1 = math.Cos(eo1)
		tem5 = 1.0 - coseo1*axnl - sineo1*aynl
		tem5 = (u - aynl*coseo1 + axnl*sineo1 - eo1) / tem5
		if math.Abs(tem5) >= 0.95 {
			if tem5 > 0.0 {
				tem5 = 0.95
			} else {
				tem5 = -0.95
			}
		}
		eo1 = eo1 + tem5
		ktr = ktr + 1
	}

	/* ------------- short period preliminary quantities ----------- */
	ecose := axnl*coseo1 + aynl*sineo1
	esine := axnl*sineo1 - aynl*coseo1
	el2 := axnl*axnl + aynl*aynl
	pl := am * (1.0 - el2)
	if pl < 0.0 {
		s.error = 4
		// sgp4fix add return
		return errors.New("pl < 0.0")
	}
	rl := am * (1.0 - ecose)
	rdotl := math.Sqrt(am) * esine / rl
	rvdotl := math.Sqrt(pl) / rl
	betal := math.Sqrt(1.0 - el2)
	temp = esine / (1.0 + betal)
	sinu := am / rl * (sineo1 - aynl - axnl*temp)
	cosu := am / rl * (coseo1 - axnl + aynl*temp)
	su := math.Atan2(sinu, cosu)
	sin2u := (cosu + cosu) * sinu
	cos2u := 1.0 - 2.0*sinu*sinu
	temp = 1.0 / pl
	temp1 := 0.5 * s.j2 * temp
	temp2 := temp1 * temp

	/* -------------- update for short period periodics ------------ */
	if s.method == 'd' {
		cosisq := cosip * cosip
		s.con41 = 3.0*cosisq - 1.0
		s.x1mth2 = 1.0 - cosisq
		s.x7thm1 = 7.0*cosisq - 1.0
	}
	mrt := rl*(1.0-1.5*temp2*betal*s.con41) + 0.5*temp1*s.x1mth2*cos2u
	su = su - 0.25*temp2*s.x7thm1*sin2u
	xnode := nodep + 1.5*temp2*cosip*sin2u
	xinc := xincp + 1.5*temp2*cosip*sinip*cos2u
	mvt := rdotl - dspaceVars.nm*temp1*s.x1mth2*sin2u/s.xke
	rvdot := rvdotl + dspaceVars.nm*temp1*(s.x1mth2*cos2u+1.5*s.con41)/s.xke

	/* --------------------- orientation vectors ------------------- */
	sinsu := math.Sin(su)
	cossu := math.Cos(su)
	snod := math.Sin(xnode)
	cnod := math.Cos(xnode)
	sini := math.Sin(xinc)
	cosi := math.Cos(xinc)
	xmx := -snod * cosi
	xmy := cnod * cosi
	ux := xmx*sinsu + cnod*cossu
	uy := xmy*sinsu + snod*cossu
	uz := sini * sinsu
	vx := xmx*cossu - cnod*sinsu
	vy := xmy*cossu - snod*sinsu
	vz := sini * cossu

	/* --------- position and velocity (in km and km/sec) ---------- */
	r[0] = (mrt * ux) * s.radiusearthkm
	r[1] = (mrt * uy) * s.radiusearthkm
	r[2] = (mrt * uz) * s.radiusearthkm
	v[0] = (mvt*ux + rvdot*vx) * vkmpersec
	v[1] = (mvt*uy + rvdot*vy) * vkmpersec
	v[2] = (mvt*uz + rvdot*vz) * vkmpersec

	// sgp4fix for decaying satellites
	if mrt < 1.0 {
		s.error = 6
		return errors.New("mrt < 1.0")
	}

	return nil
}
