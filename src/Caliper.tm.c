/*
 * This file automatically produced by /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions/mprep from:
 *	/Users/vicent/GitHubProjects/Caliper/src/Caliper.tm
 * mprep Revision 18 Copyright (c) Wolfram Research, Inc. 1990-2013
 */

#define MPREP_REVISION 18

#if CARBON_MPREP
#include <Carbon/Carbon.h>
#include <mathlink/mathlink.h>
#else
#include "mathlink.h"
#endif

int MLAbort = 0;
int MLDone  = 0;
long MLSpecialCharacter = '\0';

MLINK stdlink = 0;
MLEnvironment stdenv = 0;
#if MLINTERFACE >= 3
MLYieldFunctionObject stdyielder = (MLYieldFunctionObject)0;
MLMessageHandlerObject stdhandler = (MLMessageHandlerObject)0;
#else
MLYieldFunctionObject stdyielder = 0;
MLMessageHandlerObject stdhandler = 0;
#endif /* MLINTERFACE >= 3 */

#if DARWIN_MATHLINK && CARBON_MPREP
#define rMenuBar	1128
#define rAboutBox	1128
#define rBadSIZE	1127
#define mApple		1128
#define iAbout		1
#define mFile		1129
#define mEdit		1130

AEEventHandlerUPP handle_core_ae_upp;
ModalFilterUPP about_filter_upp;
UserItemUPP outline_hook_upp;

extern pascal OSErr handle_core_ae( const AppleEvent* event, AppleEvent* reply, long refcon);
extern int init_macintosh( void);
extern void do_about_box( void);
extern int _handle_user_event( unsigned long ticks);

pascal OSErr handle_core_ae( const AppleEvent* event, AppleEvent* reply, long refcon)
{
	DescType eventid, gottype;
	Size got;
	reply = (AppleEvent*)0; /* suppress unused warning */
	refcon = 0; /* suppress unused warning */
	if( noErr == AEGetAttributePtr(event, keyEventIDAttr, typeType, &gottype, (Ptr)&eventid, sizeof(eventid), &got)
	&& errAEDescNotFound == AEGetAttributePtr(event, keyMissedKeywordAttr, typeWildCard, &gottype, nil, 0, &got)
	){
		switch(eventid){
		case kAEQuitApplication:
			MLDone = MLAbort = 1;
		case kAEOpenApplication:
			return noErr;
		}
	}
	return errAEEventNotHandled;
}


static void set_about_item(void){
	Str255 aboutitem;
	StringHandle abouthandle;

	GetMenuItemText( GetMenuHandle(mApple), iAbout, aboutitem);
	abouthandle = NewString( aboutitem);
	if( abouthandle){
		StringPtr curApName = LMGetCurApName();
		long len = Munger( (Handle)abouthandle, 1, "MathLink\252", 9, curApName + 1, *curApName); 
		if( len > 0){
			**abouthandle = (unsigned char)len; 
			HLock( (Handle)abouthandle);
			SetMenuItemText( GetMenuHandle(mApple), iAbout, *abouthandle);
		}
		DisposeHandle( (Handle)abouthandle);
	}
}


static pascal Boolean about_filter(DialogPtr dptr, EventRecord *theEvent, short *theItem){
	if( theEvent->what == keyDown || theEvent->what == autoKey){
		unsigned char theKey = (unsigned char)(theEvent->message & charCodeMask);
		if( theKey == 0x0D || (theKey == 0x03 && !(theEvent->modifiers & controlKey))){
			short itemType;
			ControlHandle okHdl;
			Rect itemRect;
#if UNIVERSAL_INTERFACES_VERSION >= 0x0301
			unsigned long Ticks;
#else
			long Ticks;
#endif
			GetDialogItem( dptr, ok, &itemType, (Handle*) &okHdl, &itemRect);
			HiliteControl( okHdl, kControlButtonPart);
#ifdef __cplusplus
			Delay( 3, &Ticks);
#else
			Delay( 3, (void *)&Ticks);
#endif
			HiliteControl( okHdl, 0);
			*theItem = ok;
			return true;
		}
	}
	return false;
}

static pascal void outline_hook(DialogRef dptr, short theItem){
	short  itemType;
	Handle itemHdl;
	Rect itemRect;
	PenState oldpen;
	GetDialogItem( dptr, theItem, &itemType, &itemHdl, &itemRect);
	GetPenState( &oldpen);
	PenSize( 3, 3);
	FrameRoundRect( &itemRect, 16, 16);
	SetPenState( &oldpen);
}



/* edit here and in mathlink.r */
static unsigned short missing_DITL[] = { 5,
	0, 0, 76, 120, 96, 200, 0x5C, 0x30, 0x30, 0x34, 0x5C, 0x30, 0x30, 0x32, 0x4F, 0x4B, /* '\004\002', 'OK', */
	0, 0, 12, 13, 30, 28, 0x5C, 0x32, 0x31, 0x30, 0x5C, 0x30, 0x30, 0x31, 0x41, 0x5C, 0x30, /*'\210\001', 'A\0', */
	0, 0, 12, 27, 30, 96, 0x5C, 0x32, 0x31, 0x30, 0x5C, 0x30, 0x31, 0x30, 0x20, 0x4D, 0x61, 0x74, 0x68, 0x4C, 0x69, 0x6E, 0x6B, /*'\210\010', 'Ma','th','Li','nk', */
	0, 0, 12, 95, 30, 308, 0x5C, 0x32, 0x31, 0x30, 0x5C, 0x30, 0x33, 0x34, 0x5C, 0x32, 0x35, 0x32, 0x20, 0x70, 0x72, 0x6F, 0x67, 0x72, 0x61, 0x6D, 0x20, 0x67, 0x65, 0x6E, 0x65, 0x72, 0x61, 0x74, 0x65, 0x64, 0x20, 0x62, 0x79, 0x20, 0x6D, 0x70, 0x72, 0x65, 0x70, /*'\210\034', '\252','pr','og','ra','m ','ge','ne','ra','te','d ','by',' m','pr','ep', */
	0, 0, 42, 49, 56, 271, 0x5C, 0x32, 0x31, 0x30, 0x5C, 0x30, 0x35, 0x30, 0x6D, 0x70, 0x72, 0x65, 0x70, 0x5C, 0x32, 0x35, 0x31, 0x20, 0x57, 0x6F, 0x6C, 0x66, 0x72, 0x61, 0x6D, 0x20, 0x52, 0x65, 0x73, 0x65, 0x61, 0x72, 0x63, 0x68, 0x2C, 0x20, 0x49, 0x6E, 0x63, 0x2E, 0x20, 0x31, 0x39, 0x39, 0x30, 0x2D, 0x32, 0x30, 0x30, 0x32, /*'\210\050', 'mp','re','p ','\251 ','Wo','lf','ra','m','Re','se','ar','ch',', ','In','c.',' 1','99','0-','20','02', */ /* 1990-2002 */
	0, 0, 170, 10, 190, 30, 0x5C, 0x32, 0x30, 0x30, 0x5C, 0x30, 0x30, 0x30 /*'\200\000' */
};


int init_macintosh( void)
{
	static int initdone = 0;
	Handle menuBar;
	long attributes;

	/* semaphore required for preemptive threads */
	if( initdone) return initdone == 1;
	initdone = -1;

	/* should I check for MLNK resource too as launch-filtering is done based on it? */
	/* too late--since I'm running there likely wasn't a problem (this time anyway). */
	
	menuBar = GetNewMBar(rMenuBar);
	if( menuBar){
		SetMenuBar(menuBar);
		DisposeHandle(menuBar);
	}else{
		MenuHandle am, fm, em;
		am = NewMenu( mApple, (unsigned char*)"\001\024");
		fm = NewMenu( mFile, (unsigned char*)"\004File");
		em = NewMenu( mEdit, (unsigned char*)"\004Edit");
		if( !am || !fm || !em) return 0;
		AppendMenu( am, (unsigned char*)"\022About MathLink\252\311;-");
                DisableMenuItem(am, 0);
		InsertMenu( am, 0);
		AppendMenu( fm, (unsigned char*)"\006Quit/Q");
		InsertMenu( fm, 0);
		AppendMenu( em, (unsigned char*)"\043Undo/Z;-;Cut/X;Copy/C;Paste/V;Clear");
                DisableMenuItem(em, 0);
		InsertMenu( em, 0);
	}

	AppendResMenu( GetMenuHandle(mApple), 'DRVR');
	set_about_item();
	DrawMenuBar();
	about_filter_upp =  NewModalFilterUPP( about_filter);
	outline_hook_upp = NewUserItemUPP( outline_hook);
	if( Gestalt( gestaltAppleEventsAttr, &attributes) == noErr
	&& ((1 << gestaltAppleEventsPresent) & attributes)){
		handle_core_ae_upp = NewAEEventHandlerUPP( handle_core_ae);
		(void) AEInstallEventHandler( kCoreEventClass, typeWildCard, handle_core_ae_upp, 0, false);
	}else{
		return 0; /* this may be too strong since I am, after all, running. */
	}

	initdone = 1;
	return initdone == 1;
}

void do_about_box( void)
{
	GrafPtr oldPort;
	DialogPtr dptr;
	short item, itemType;
	Handle itemHdl;
	Rect itemRect;

	dptr = GetNewDialog( rAboutBox, nil, (WindowPtr)-1L);
	
	if( dptr == (DialogPtr)0){
		Handle items = NewHandle( sizeof(missing_DITL));
		static Rect bounds = {40, 20, 150, 340};

		if( ! items) return;
		BlockMove( missing_DITL, *items, sizeof(missing_DITL));

		dptr = NewColorDialog( nil, &bounds, (unsigned char*)"\005About",
					false, dBoxProc, (WindowPtr)-1L, false, 0, items);
                }
	
	if( dptr == (DialogPtr)0) return;
	GetPort (&oldPort);
	SetPort (GetDialogPort(dptr));
	GetDialogItem( dptr, ok, &itemType, &itemHdl, &itemRect);
	InsetRect( &itemRect, -4, -4);
	SetDialogItem( dptr, 6, userItem + itemDisable, (Handle)outline_hook_upp, &itemRect);

	FlushEvents( everyEvent, 0);
        ShowWindow( GetDialogWindow(dptr));

	do {
		ModalDialog( about_filter_upp, &item);
	} while( item != ok);

	DisposeDialog( dptr);
	SetPort( oldPort);
}

int _handle_user_event( unsigned long ticks)
{
	EventRecord event;

	if( WaitNextEvent(everyEvent, &event, ticks, nil)){
		long      menuResult = 0;
		short     menuID, menuItem;
		WindowPtr window;
		Str255    daName;

		switch ( event.what ) {
		case mouseDown:
			if( FindWindow(event.where, &window) == inMenuBar)
				menuResult = MenuSelect(event.where);
			break;
		case keyDown:
			if( event.modifiers & cmdKey )
				menuResult = MenuKey((short)event.message & charCodeMask);
			break;
		case kHighLevelEvent:
			AEProcessAppleEvent(&event);
			break;
		}

		menuID = HiWord(menuResult);
		menuItem = LoWord(menuResult);
		switch ( menuID ) {
		case mFile:
			MLDone = MLAbort = 1;
			break;
		case mApple:
			switch ( menuItem ) {
			case iAbout:
				do_about_box();
				break;
			default:
				GetMenuItemText(GetMenuHandle(mApple), menuItem, daName);
				break;
			}
			HiliteMenu(0);
		}
	}
	return MLDone;
}

#if MLINTERFACE >= 3
MLYDEFN( int, MLDefaultYielder, ( MLINK mlp, MLYieldParameters yp))
#else
MLYDEFN( devyield_result, MLDefaultYielder, ( MLINK mlp, MLYieldParameters yp))
#endif /* MLINTERFACE >= 3 */
{
	mlp = mlp; /* suppress unused warning */
#if MLINTERFACE >= 3
	return (int)_handle_user_event( MLSleepYP(yp));
#else
	return _handle_user_event( MLSleepYP(yp));
#endif /* MLINTERFACE >= 3 */
}

#endif /* (DARWIN_MATHLINK && CARBON_MPREP */

/********************************* end header *********************************/


# line 1576 "/Users/vicent/GitHubProjects/Caliper/src/Caliper.tm"
#include "mathlink.h"
#include "ftypes.h"
#include <stdio.h>
#include <unistd.h>

extern double f90ewfactors_(int* nf, double* Q, double* Mz, double* GammaZ, double* sin2ThetaW, double* res);

static void ewfactors(int nf, double Q, double Mz, double GammaZ, double sin2ThetaW){
  double res[2];

   f90ewfactors_(&nf, &Q, &Mz, &GammaZ, &sin2ThetaW, res);

   MLPutRealList(stdlink, res, 2);
   MLEndPacket(stdlink);
}

extern double f90kernels_(int* n, double* width, double* w, double* mu, double* p, double* result);

static void kernels(int n, double width, double w, double mu, double p){
  double res[n+1];

   f90kernels_(&n, &width, &w, &mu, &p, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);
}

extern double f90gammaderlist_(int* n, double* w, double* result);

static void gammaderlist(int n, double w){
  double res[n+1];

   f90gammaderlist_(&n, &w, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);
}

extern double f90polygamma_(int* n, double* w, double* result);

static void polygamma(int n, double w){
  double res[n+1];

   f90polygamma_(&n, &w, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);
}

extern double f90nglkernels_(int* n, int* n1, int*n2, double* width, double* w, double* mu, double* p, double* result);

static void nglkernels(int n, int n1, int n2, double width, double w[], long wlen, double mu, double p[], long plen){
  double res[n];

   f90nglkernels_(&n, &n1, &n2, &width, w, &mu, p, res);

   MLPutRealList(stdlink, res, n);

   MLEndPacket(stdlink);
}

extern double f90nglintegral_(int* nf, int* pow, double* w1, double* w2, double* res);

static double nglintegral(int nf, int pow, double w1, double w2){
  double res;

   f90nglintegral_(&nf, &pow, &w1, &w2, &res);

   return res;
}

extern double f90ngldoubleintegral_(int* nf, int* pow, double* w1, double* w2, double* r, double* res);

static double ngldoubleintegral(int nf, int pow, double w1, double w2, double r){
  double res;

   f90ngldoubleintegral_(&nf, &pow, &w1, &w2, &r, &res);

   return res;
}

extern double f90singularhjmpiece_(char const* hard, char const* gap, char const* space, char const* cum,
 int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf, double* j3,
 double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, int* c,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* res);

static double singularhjmpiece(char const* hard, char const* gap, char const* space, char const* cum,
 int orderAlpha, int runAlpha, int order, int run, int isoft, int nf, double j3, double s3,
 double s31, double s32, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
 double muJ, double muS, double R, double mu, int c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau){
  double res;

f90singularhjmpiece_(hard, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q,
 &muH, &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singulardoublepiece_(char const* hard, char const* gap, char const* space, char const* cum1,
 char const* cum2, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft,
 int* nf, double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ,
 double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 int* c, double* lambda, double* R0, double* mu0, double* delta0, double* h, double* rho1,
 double* rho2, double* res);

static double singulardoublepiece(char const* hard, char const* gap, char const* space, char const* cum1,
 char const* cum2, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, int c[], long len, double lambda,
 double R0, double mu0, double delta0, double h, double rho1, double rho2){
  double res;

f90singulardoublepiece_(hard, gap, space, cum1, cum2, &orderAlpha, &runAlpha, &order, &run,
 &isoft, &nf, &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
 &muLambda, &Q, &muH, &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &rho1,
 &rho2, &res);

return res;

}

extern double f90legendrelist_(int* n, int* k, double* x, double* res);

static void legendrelist(int n, int k, double x){
  double res[n + 1];

   f90legendrelist_(&n, &k, &x, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);

}

extern double f90qlegendrelist_(int* n, double* x, double* res);

static void qlegendrelist(int n, double x){
  double res[n + 1];

   f90qlegendrelist_(&n, &x, res);

   MLPutRealList(stdlink, res, n + 1);
   MLEndPacket(stdlink);

}

extern double f90singularhjm_(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf,
 double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static double singularhjm(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, double c[], long len, int clen,
 double lambda, double R0, double mu0, double delta0, double h, double tau){
  double res;

f90singularhjm_(hard, setup, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ,
 &muS, &R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singularhjm1d_(char const* hard, char const* gap,
 char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf,
 double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static double singularhjm1d(char const* hard, char const* gap,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, double c[], long len,
 double lambda, double R0, double mu0, double delta0, double h, double tau){
  double res;

 int clen = len;

f90singularhjm1d_(hard, gap, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q,
 &muH, &muJ, &muS, &R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singularhjm1dpiece_(char const* hard, char const* gap,
 char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run, int* isoft, int* nf,
 double* j3, double* s3, double* s31, double* s32, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda, double* Q, double* muH,
 double* muJ, double* muS, double* R, double* mu, int* c, double* lambda,
 double* R0, double* mu0, double* delta0, double* h, double* tau, double* res);

static double singularhjm1dpiece(char const* hard, char const* gap,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int isoft, int nf,
 double j3, double s3, double s31, double s32, double G3, double mZ, double aMz, double mT,
 double muT, double mB, double muB, double mC, double muC, double muLambda, double Q,
 double muH, double muJ, double muS, double R, double mu, int c[], long len, double lambda,
 double R0, double mu0, double delta0, double h, double tau){
  double res;

f90singularhjm1dpiece_(hard, gap, cum, &orderAlpha, &runAlpha, &order, &run, &isoft, &nf,
 &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q,
 &muH, &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singulardouble_(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum1, char const* cum2, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* isoft, int* nf, double* j3, double* s3, double* s31, double* s32, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* mu, double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0,
 double* h, double* rho1, double* rho2, double* res);

static double singulardouble(char const* hard, char const* setup, char const* gap, char const* space,
 char const* cum1, char const* cum2, int orderAlpha, int runAlpha, int order, int run,
 int isoft, int nf, double j3, double s3, double s31, double s32,double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double muH, double muJ, double muS, double R, double mu,
 double c[], long len, int clen, double lambda, double R0, double mu0, double delta0,
 double h, double rho1, double rho2){
  double res;

f90singulardouble_(hard, setup, gap, space, cum1, cum2, &orderAlpha, &runAlpha, &order, &run,
 &isoft, &nf, &j3, &s3, &s31, &s32, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
 &muLambda, &Q, &muH, &muJ, &muS, &R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h,
 &rho1, &rho2, &res);

return res;

}

extern double f90profiles_(double* Q, double* mu0, double* R0, double* n0, double* n1,
double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH, double* eS,
double* eJ, double* eR, int* ns, double* tau, double* res);

static void profiles(double Q, double mu0, double R0, double n0, double n1, double t2,
double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
double eR, int ns, double tau){
  double result[5];

   f90profiles_(&Q, &mu0, &R0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR,
   &ns, &tau, result);

   MLPutRealList(stdlink, result, 5);

   MLEndPacket(stdlink);
}

extern double f90profilesmass_(double* Q, double* beta, double* mu0, double* delLamb,
double* R0, double* n0, double* delta0, double* n1, double* delta1, double* t2, double* ts,
double* slope, double* cnt, double* eH, double* eS, double* eJ, double* mass, double* muM,
int* ns, char const* def, char const* EShape, double* tau, double* res);

static void profilesmass(double Q, double beta, double mu0, double delLamb, double R0,
double n0, double delta0, double n1, double delta1, double t2, double ts, double slope,
double cnt, double eH, double eS, double eJ, double mass, double muM, int ns,
char const* def, char const* EShape, double tau){
  double result[6];

   f90profilesmass_(&Q, &beta, &mu0, &delLamb, &R0, &n0, &delta0, &n1, &delta1, &t2,
   &ts, &slope, &cnt, &eH, &eS, &eJ, &mass, &muM, &ns, def, EShape, &tau, result);

   MLPutRealList(stdlink, result, 6);

   MLEndPacket(stdlink);
}

extern double f90massintersection_(double* Q, double* beta, double* mu0, double* delLamb,
double* n0, double* delta0, double* n1, double* delta1, double* t2, double* ts,
double* slope, double* cnt, double* eH, double* eS, double* eJ, double* mass, double* muM,
char const* def, char const* EShape, double* res);

static void massintersection(double Q, double beta, double mu0, double delLamb, double n0,
double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
double eH, double eS, double eJ, double mass, double muM, char const* def,
char const* EShape){
  double result[2];

   f90massintersection_(&Q, &beta, &mu0, &delLamb, &n0, &delta0, &n1, &delta1, &t2,
   &ts, &slope, &cnt, &eH, &eS, &eJ, &mass, &muM, def, EShape, result);

   MLPutRealList(stdlink, result, 2);

   MLEndPacket(stdlink);
}

extern double f90gammar_(char const* str, int* nf, double* result);

static void gammar(char const* str, int nf){
  double result[3];

   f90gammar_(str, &nf, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90masslessprof_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* tau, double* res);

static double masslessprof(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double tau){
  double res;
  int clen = len;

f90masslessprof_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &res);

return res;

}

extern double f90findorigin_(char const* shape, char const* gap, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* mu0, double* Rat0, double* n0, double* n1, double* t2, double* tR,
 double* ts, double* slope, double* cnt, double* eH, double* eS, double* eR,
 double* R0, double* muR0, double* delta0, double* h, double* res);

static double findorigin(char const* shape, char const* gap, int orderAlpha, int runAlpha,
 int order, int run, int nf, double mZ, double aMz, double mT, double muT, double mB,
 double muB, double mC, double muC, double muLambda, double Q, double mu0, double Rat0,
 double n0, double n1, double t2, double tR, double ts, double slope, double cnt,
 double eH, double eS, double eR, double R0, double muR0, double delta0, double h){
  double res;

f90findorigin_(shape, gap, &orderAlpha, &runAlpha, &order, &run, &nf, &mZ, &aMz, &mT,
&muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope,
&cnt, &eH, &eS, &eR, &R0, &muR0, &delta0, &h, &res);

return res;

}

extern double f90masslessprofpiece_(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* tau, double* res);

static void masslessprofpiece(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double tau){
  double res[(clen + 1) * (clen + 2)/2];

f90masslessprofpiece_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslessprofpiecelist_(char const* terms, char const* hard,
 char const* shape, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* taulist, int* taulen, double* res);

static void masslessprofpiecelist(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double taulist[], long tlen){
  int taulen = tlen;
  double res[tlen * (clen + 1) * (clen + 2)/2];

f90masslessprofpiecelist_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, tlen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslesspiecebin_(char const* terms, char const* hard,
 char const* shape, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* taulist, int* taulen, double* res);

static void masslesspiecebin(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double taulist[], long tlen, int taulen){
  double res[taulen * (clen + 1) * (clen + 2)/2];

f90masslesspiecebin_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, taulen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslessprofdiffpiece_(char const* terms, char const* hard,
char const* shape, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* tau, double* tau2, double* res);

static void masslessprofdiffpiece(char const* terms, char const* hard, char const* shape,
 char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, int clen, double lambda, double R0, double muR0, double delta0,
 double h, double tau, double tau2){
  double res[(clen + 1) * (clen + 2)/2];

f90masslessprofdiffpiece_(terms, hard, shape, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &tau2, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90masslessproflist_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* taulist, int* tlen,
 double* result);

static void masslessproflist(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double taulist[], long taulen){
  double result[taulen];
  int clen = len;
  int tlen = taulen;

f90masslessproflist_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &tlen, result);

   MLPutRealList(stdlink, result, taulen);

   MLEndPacket(stdlink);
}

extern double f90masslessbinlist_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* taulist, int* tlen,
 double* result);

static void masslessbinlist(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double taulist[], long taulen, int tlen){
  double result[tlen];
  int clen = len;

f90masslessbinlist_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, taulist, &tlen, result);

   MLPutRealList(stdlink, result, tlen);

   MLEndPacket(stdlink);
}

extern double f90masslessprofdiff_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int* orderAlpha,
 int* runAlpha, int* order, int* run, int* nf, double* j3, double* s3, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda, double* Q, double* mu0, double* Rat0, double* n0,
 double* n1, double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH,
 double* eS, double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* delta0, double* h, double* tau, double* tau2,
 double* res);

static double masslessprofdiff(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, char const* cum, int orderAlpha,
 int runAlpha, int order, int run, int nf, double j3, double s3, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda, double Q, double mu0, double Rat0, double n0, double n1, double t2,
 double tR, double ts, double slope, double cnt, double eH, double eS, double eJ,
 double eR, int ns, double c[], long len, double lambda, double R0, double muR0,
 double delta0, double h, double tau, double tau2){
  double res;
  int clen = len;

f90masslessprofdiff_(terms, hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90masslessmoment_(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, int* orderAlpha, int* runAlpha,
 int* order, int* run, int* nf, double* j3, double* s3, double* G3, double* mZ,
 double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda, double* Q, double* mu0, double* Rat0, double* n0, double* n1,
 double* t2, double* tR, double* ts, double* slope, double* cnt, double* eH, double* eS,
 double* eJ, double* eR, int* ns, double* c, int* clen, double* lambda, double* R0,
 double* muR0, double* delta0, double* h, double* tau, double* tau2, int* pow,
 double* res);

static double masslessmoment(char const* terms, char const* hard, char const* shape,
 char const* setup, char const* gap, char const* space, int orderAlpha, int runAlpha,
 int order, int run, int nf, double j3, double s3, double G3, double mZ, double aMz,
 double mT, double muT, double mB, double muB, double mC, double muC, double muLambda,
 double Q, double mu0, double Rat0, double n0, double n1, double t2, double tR, double ts,
 double slope, double cnt, double eH, double eS, double eJ, double eR, int ns, double c[],
 long len, double lambda, double R0, double muR0, double delta0, double h, double tau,
 double tau2, int pow){
  double res;
  int clen = len;

f90masslessmoment_(terms, hard, shape, setup, gap, space, &orderAlpha, &runAlpha,
&order, &run, &nf, &j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda,
&Q, &mu0, &Rat0, &n0, &n1, &t2, &tR, &ts, &slope, &cnt, &eH, &eS, &eJ, &eR, &ns, c, &clen,
&lambda, &R0, &muR0, &delta0, &h, &tau, &tau2, &pow, &res);

return res;

}

extern double f90singular_(char const* hard, char const* shape, char const* setup, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static double singular(char const* hard, char const* shape, char const* setup, char const* gap,
char const* space, char const* cum, int orderAlpha, int runAlpha, int order, int run,
int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
double muJ, double muS, double R, double mu, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double tau){
  double res;
  int clen = len;

f90singular_(hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ, &muS,
&R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singularlist_(char const* hard, char const* shape, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau, double* res);

static void singularlist(char const* hard, char const* shape, char const* gap,
char const* space, char const* cum, int orderAlpha, int runAlpha, int order, int run,
int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
double muJ, double muS, double R, double mu, int clen, double lambda,
double R0, double mu0, double delta0, double h, double tau){
  double res[(clen + 1) * (clen + 2)/2];

f90singularlist_(hard, shape, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ, &muS,
&R, &mu, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);

}


extern double f90singulardiff_(char const* hard, char const* shape, char const* setup, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* tau1, double* tau2, double* res);

static double singulardiff(char const* hard, char const* shape, char const* setup, char const* gap,
char const* space, char const* cum, int orderAlpha, int runAlpha, int order, int run,
int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
double mB, double muB, double mC, double muC, double muLambda, double Q, double muH,
double muJ, double muS, double R, double mu, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double tau1, double tau2){
  double res;
  int clen = len;

f90singulardiff_(hard, shape, setup, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ, &muS,
&R, &mu, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau1, &tau2, &res);

return res;

}

extern double f90massiveprof_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* res);

static double massiveprof(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau){
  double res;
  int clen = len;

f90massiveprof_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &res);

return res;

}

extern double f90massorigin_(char const* shape, char const* Eshape, char const* gap,
 char const* scheme, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass,
 int* order, int* run, int* nf, double* mZ, double* aMz, double* mT, double* muT,
 double* mB, double* muB, double* mC, double* muC, double* muLambda1, double* muLambda2,
 double* Q, double* beta, double* mu0, double* deltaLambda, double* Rat0, double* n0,
 double* delta0, double* n1, double* delta1, double* t2, double* ts, double* slope,
 double* cnt, double* eH, double* eS, double* eJ, double* mass, double* muM,
 double* R0, double* muR0, double* del0, double* h, double* res);

static double massorigin(char const* shape, char const* Eshape, char const* gap,
 char const* scheme, int orderAlpha, int runAlpha, int orderMass, int runMass, int order,
 int run, int nf, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda1, double muLambda2, double Q, double beta,
 double mu0, double deltaLambda, double Rat0, double n0, double delta0, double n1,
 double delta1, double t2, double ts, double slope, double cnt, double eH, double eS,
 double eJ, double mass, double muM, double R0, double muR0, double del0, double h){
  double res;

f90massorigin_(shape, Eshape, gap, scheme, &orderAlpha, &runAlpha, &orderMass, &runMass,
 &order, &run, &nf, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2,
 &Q, &beta, &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt,
 &eH, &eS, &eJ, &mass, &muM, &R0, &muR0, &del0, &h, &res);

return res;

}

extern double f90massiveprofpiece_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* res);

static void massiveprofpiece(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau){
  double res[(clen + 1) * (clen + 2)/2];

f90massiveprofpiece_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massiveprofpiecelist_(char const* terms, char const* hard,
 char const* shape, char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, int* taulen, double* res);

static void massiveprofpiecelist(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau[], long tlen){
  int taulen = tlen;
  double res[tlen * (clen + 1) * (clen + 2)/2];

f90massiveprofpiecelist_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tau, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, tlen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massivepiecebin_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, int* taulen, double* res);

static void massivepiecebin(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau[], long tlen, int taulen){
  double res[taulen * (clen + 1) * (clen + 2)/2];

f90massivepiecebin_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tau, &taulen, res);

   MLPutFunction(stdlink, "Partition", 2 );
   MLPutRealList(stdlink, res, taulen * (clen + 1) * (clen + 2)/2);
   MLPutInteger(stdlink, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massiveprofdiffpiece_(char const* terms, char const* hard,
 char const* shape, char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* tau2, double* res);

static void massiveprofdiffpiece(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double xi, double xiB,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 int clen, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau, double tau2){
  double res[(clen + 1) * (clen + 2)/2];

f90massiveprofdiffpiece_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &tau2, res);

   MLPutRealList(stdlink, res, (clen + 1) * (clen + 2)/2);

   MLEndPacket(stdlink);
}

extern double f90massiveprofdiff_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* tau2, double* res);

static double massiveprofdiff(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau, double tau2){
  double res;
  int clen = len;

f90massiveprofdiff_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &tau2, &res);

return res;

}

extern double f90massiveproflist_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tauList, int* taulen, double* res);

static void massiveproflist(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tauList[], long ctau){
  int clen = len;
  int taulen = ctau;
  double result[taulen];

f90massiveproflist_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tauList, &taulen, result);

   MLPutRealList(stdlink, result, taulen);

   MLEndPacket(stdlink);

}

extern double f90massivebinlist_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tauList, int* taulen, double* res);

static void massivebinlist(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tauList[], long ctau, int taulen){
  int clen = len;
  double result[taulen];

f90massivebinlist_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,
 &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3,
 &G3,  &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 tauList, &taulen, result);

   MLPutRealList(stdlink, result, taulen);

   MLEndPacket(stdlink);

}

extern double f90massivediffprof_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* tau2, double* res);

static double massivediffprof(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* cum, char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau, double tau2){
  double res;
  int clen = len;

f90massivediffprof_(terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs,
 current, &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3,
 &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &tau2, &res);

return res;

}

extern double f90massivemoment_(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* scheme, char const* abs, char const* current, double* xi,
 double* xiB, int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order,
 int* run, int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* beta, double* mu0,
 double* deltaLambda, double* Rat0, double* n0, double* delta0, double* n1, double* delta1,
 double* t2, double* ts, double* slope, double* cnt, double* eH, double* eS, double* eJ,
 double* mass, double* muM, int* ns, double* width, double* c, int* clen, double* lambda,
 double* R0, double* muR0, double* del0, double* h, double* gammaZ, double* sinW,
 double* tau, double* tau2, int* pow, double* res);

static double massivemoment(char const* terms, char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space,
 char const* scheme, char const* abs, char const* current, double xi,
 double xiB, int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run,
 int nf, double j3, double s3, double G3, double mZ, double aMz, double mT, double muT,
 double mB, double muB, double mC, double muC, double muLambda1, double muLambda2,
 double Q, double beta, double mu0, double deltaLambda, double Rat0, double n0,
 double delta0, double n1, double delta1, double t2, double ts, double slope, double cnt,
 double eH, double eS, double eJ, double mass, double muM, int ns, double width,
 double c[], long len, double lambda, double R0, double muR0, double del0, double h,
 double gammaZ, double sinW, double tau, double tau2, int pow){
  double res;
  int clen = len;

f90massivemoment_(terms, hard, shape, Eshape, setup, gap, space, scheme, abs,
 current, &xi, &xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3,
 &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &beta,
 &mu0, &deltaLambda, &Rat0, &n0, &delta0, &n1, &delta1, &t2, &ts, &slope, &cnt, &eH, &eS,
 &eJ, &mass, &muM, &ns, &width, c, &clen, &lambda, &R0, &muR0, &del0, &h, &gammaZ, &sinW,
 &tau, &tau2, &pow, &res);

return res;

}

extern double f90singularmass_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* res);

static double singularmass(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
char const* abs, char const* current, double xi, double xiB, int orderAlpha, int runAlpha,
int orderMass, int runMass, int order, int run, int nf, double j3, double s3, double G3,
double mZ, double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;
  int clen = len;

f90singularmass_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW,
&tau, &res);

return res;

}

extern double f90singularmassdiff_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* tau2, double* res);

static double singularmassdiff(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
char const* abs, char const* current, double xi, double xiB, int orderAlpha, int runAlpha,
int orderMass, int runMass, int order, int run, int nf, double j3, double s3, double G3,
double mZ, double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau,
double tau2){
  double res;
  int clen = len;

f90singularmassdiff_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &clen, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW,
&tau, &tau2, &res);

return res;

}

extern double f90massnondist_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf,
 double* G3, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH,
 double* muJ, double* muS, double* R, double* Rmass, double* muM, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* res);

static double massnondist(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run, int nf,
 double G3, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda1, double muLambda2, double Q, double muH, double muJ,
 double muS, double R, double Rmass, double muM, double mu, double c[], long len,
 double lambda, double R0, double mu0, double delta0, double h, double tau){
  double res;
  int clen = len;

f90massnondist_(hard, shape, Eshape, setup, gap, space, cum, scheme, &orderAlpha,
 &runAlpha, &orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB,
 &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu,
 c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90massnondistdiff_(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf,
 double* G3, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH,
 double* muJ, double* muS, double* R, double* Rmass, double* muM, double* mu, double* c,
 int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* tau2, double* res);

static double massnondistdiff(char const* hard, char const* shape, char const* Eshape,
 char const* setup, char const* gap, char const* space, char const* cum, char const* scheme,
 int orderAlpha, int runAlpha, int orderMass, int runMass, int order, int run, int nf,
 double G3, double mZ, double aMz, double mT, double muT, double mB, double muB, double mC,
 double muC, double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS,
 double R, double Rmass, double muM, double mu, double c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau, double tau2){
  double res;
  int clen = len;

f90massnondistdiff_(hard, shape, Eshape, setup, gap, space, cum, scheme, &orderAlpha,
&runAlpha, &orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB,
&muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu,
c, &clen, &lambda, &R0, &mu0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90massnondistpiece_(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int* orderAlpha,
 int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH, double* muJ,
 double* muS, double* R, double* Rmass, double* muM, double* mu, int* c, double* lambda,
 double* R0, double* mu0, double* delta0, double* h, double* tau, double* res);

static double massnondistpiece(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int orderAlpha,
 int runAlpha, int orderMass, int runMass, int order, int run, int nf, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
 double Rmass, double muM, double mu, int c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau){
  double res;

f90massnondistpiece_(hard, shape, Eshape, gap, space, cum, scheme, &orderAlpha, &runAlpha,
 &orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
 &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu, c, &lambda, &R0,
 &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90massnondistdiffpiece_(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int* orderAlpha,
 int* runAlpha, int* orderMass, int* runMass, int* order, int* run, int* nf, double* G3,
 double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
 double* muC, double* muLambda1, double* muLambda2, double* Q, double* muH, double* muJ,
 double* muS, double* R, double* Rmass, double* muM, double* mu, int* c, double* lambda,
 double* R0, double* mu0, double* delta0, double* h, double* tau, double* tau2,
 double* res);

static double massnondistdiffpiece(char const* hard, char const* shape, char const* Eshape,
 char const* gap, char const* space, char const* cum, char const* scheme, int orderAlpha,
 int runAlpha, int orderMass, int runMass, int order, int run, int nf, double G3, double mZ,
 double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
 double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
 double Rmass, double muM, double mu, int c[], long len, double lambda, double R0,
 double mu0, double delta0, double h, double tau, double tau2){
  double res;

f90massnondistdiffpiece_(hard, shape, Eshape, gap, space, cum, scheme, &orderAlpha, &runAlpha,
&orderMass, &runMass, &order, &run, &nf, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
&muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R, &Rmass, &muM, &mu, c, &lambda,
&R0, &mu0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90singularmasspiece_(char const* hard, char const* shape,
 char const* Eshape, char const* setup, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* res);

static double singularmasspiece(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme, char const* abs,
char const* current, double xi, double xiB, int orderAlpha, int runAlpha, int orderMass,
int runMass, int order, int run, int nf, double j3, double s3, double G3, double mZ,
double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, int c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;

f90singularmasspiece_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90singularmassdiffpiece_(char const* hard, char const* shape,
 char const* setup, char const* Eshape, char const* gap, char const* space, char const* cum,
 char const* scheme, char const* abs, char const* current, double* xi, double* xiB,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda1,
 double* muLambda2, double* Q, double* muH, double* muJ, double* muS, double* R,
 double* Rmass, double* muM, double* mu, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* tau2, double* res);

static double singularmassdiffpiece(char const* hard, char const* shape, char const* Eshape,
char const* setup, char const* gap, char const* space, char const* cum, char const* scheme, char const* abs,
char const* current, double xi, double xiB, int orderAlpha, int runAlpha, int orderMass,
int runMass, int order, int run, int nf, double j3, double s3, double G3, double mZ,
double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double muLambda1, double muLambda2, double Q, double muH, double muJ, double muS, double R,
double Rmass, double muM, double mu, double width, int c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau,
double tau2){
  double res;

f90singularmassdiffpiece_(hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, &xi,
&xiB, &orderAlpha, &runAlpha, &orderMass, &runMass, &order, &run, &nf, &j3, &s3, &G3, &mZ,
&aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda1, &muLambda2, &Q, &muH, &muJ, &muS, &R,
&Rmass, &muM, &mu, &width, c, &lambda, &R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2,
&res);

return res;

}

extern double f90nsmass_(char const* shape, char const* setup, char const* gap,
 char const* cum, char const* scheme, char const* abs, char const* current,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muBottom,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* mu,
 double* muM, double* muB, double* muS, double* R, double* Rmass, double* width,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* gammaZ, double* sinW, double* tau, double* res);

static double nsmass(char const* shape, char const* setup, char const* gap,
char const* cum, char const* scheme, char const* abs, char const* current, int orderAlpha,
int runAlpha, int order, int run, int orderMass, int runMass, int nf, double mZ,
double aMz, double mT, double muT, double mB, double muBottom, double mC, double muC,
double muLambda1, double muLambda2, double Q, double mu, double muM, double muB,
double muS, double R, double Rmass, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;
  int clen = len;

f90nsmass_(shape, setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90nsmassdiff_(char const* shape, char const* setup, char const* gap,
 char const* cum, char const* scheme, char const* abs, char const* current,
 int* orderAlpha, int* runAlpha, int* orderMass, int* runMass, int* order, int* run,
 int* nf, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muBottom,
 double* mC, double* muC, double* muLambda1, double* muLambda2, double* Q, double* mu,
 double* muM, double* muB, double* muS, double* R, double* Rmass, double* width,
 double* c, int* clen, double* lambda, double* R0, double* mu0, double* delta0, double* h,
 double* gammaZ, double* sinW, double* tau, double* tau2, double* res);

static double nsmassdiff(char const* shape, char const* setup, char const* gap,
char const* cum, char const* scheme, char const* abs, char const* current, int orderAlpha,
int runAlpha, int order, int run, int orderMass, int runMass, int nf, double mZ,
double aMz, double mT, double muT, double mB, double muBottom, double mC, double muC,
double muLambda1, double muLambda2, double Q, double mu, double muM, double muB,
double muS, double R, double Rmass, double width, double c[], long len, double lambda,
double R0, double mu0, double delta0, double h, double gammaZ, double sinW, double tau,
double tau2){
  double res;
  int clen = len;

f90nsmassdiff_(shape, setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2, &res);

return res;

}

extern double f90hjmnsmass_(char const* setup, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* res);

static double hjmnsmass(char const* setup, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, double c[], long len, int clen, double lambda, double R0,
double mu0, double delta0, double h, double gammaZ, double sinW, double tau){
  double res;

f90hjmnsmass_(setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90hjmnsmassdiff_(char const* setup, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, double* c, int* clen,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* gammaZ,
 double* sinW, double* tau, double* tau2, double* res);

static double hjmnsmassdiff(char const* setup, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, double c[], long len, int clen, double lambda, double R0,
double mu0, double delta0, double h, double gammaZ, double sinW, double tau, double tau2){
  double res;

f90hjmnsmassdiff_(setup, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda2, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &clen, &lambda,
&R0, &mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2, &res);

return res;

}

extern double f90nsmasspiece_(char const* shape, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* res);

static double nsmasspiece(char const* shape, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, int c[], long len, double lambda, double R0, double mu0,
double delta0, double h, double gammaZ, double sinW, double tau){
  double res;

f90nsmasspiece_(shape, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda1, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &lambda, &R0,
&mu0, &delta0, &h, &gammaZ, &sinW, &tau, &res);

return res;

}

extern double f90nsmassdiffpiece_(char const* shape, char const* gap, char const* cum,
 char const* scheme, char const* abs, char const* current, int* orderAlpha, int* runAlpha,
 int* orderMass, int* runMass, int* order, int* run, int* nf, double* mZ, double* aMz,
 double* mT, double* muT, double* mB, double* muBottom, double* mC, double* muC,
 double* muLambda1, double* muLambda2, double* Q, double* mu, double* muM, double* muB,
 double* muS, double* R, double* Rmass, double* width, int* c, double* lambda, double* R0,
 double* mu0, double* delta0, double* h, double* gammaZ, double* sinW, double* tau,
 double* tau2, double* res);

static double nsmassdiffpiece(char const* shape, char const* gap, char const* cum,
char const* scheme, char const* abs, char const* current, int orderAlpha, int runAlpha,
int order, int run, int orderMass, int runMass, int nf, double mZ, double aMz, double mT,
double muT, double mB, double muBottom, double mC, double muC, double muLambda1,
double muLambda2, double Q, double mu, double muM, double muB, double muS, double R,
double Rmass, double width, int c[], long len, double lambda, double R0, double mu0,
double delta0, double h, double gammaZ, double sinW, double tau, double tau2){
  double res;

f90nsmassdiffpiece_(shape, gap, cum, scheme, abs, current, &orderAlpha, &runAlpha, &order,
&run, &orderMass, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB, &muBottom, &mC, &muC,
&muLambda1, &muLambda1, &Q, &mu, &muM, &muB, &muS, &R, &Rmass, &width, c, &lambda, &R0,
&mu0, &delta0, &h, &gammaZ, &sinW, &tau, &tau2, &res);

return res;

}

extern double f90singularpiece_(char const* hard, char const* shape, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, int* c,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* res);

static double singularpiece(char const* hard, char const* shape, char const* gap, char const* space,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int nf, double j3,
 double s3, double G3, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda, double Q, double muH, double muJ, double muS,
 double R, double mu, int c[], long len, double lambda, double R0, double mu0, double delta0,
 double h, double tau){
  double res;

f90singularpiece_(hard, shape, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf, &j3,
             &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH,
             &muJ, &muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &res);

return res;

}

extern double f90singulardiffpiece_(char const* hard, char const* shape, char const* gap,
 char const* space, char const* cum, int* orderAlpha, int* runAlpha, int* order, int* run,
 int* nf, double* j3, double* s3, double* G3, double* mZ, double* aMz, double* mT,
 double* muT, double* mB, double* muB, double* mC, double* muC, double* muLambda,
 double* Q, double* muH, double* muJ, double* muS, double* R, double* mu, int* c,
 double* lambda, double* R0, double* mu0, double* delta0, double* h, double* tau,
 double* tau2, double* res);

static double singulardiffpiece(char const* hard, char const* shape, char const* gap, char const* space,
 char const* cum, int orderAlpha, int runAlpha, int order, int run, int nf, double j3,
 double s3, double G3, double mZ, double aMz, double mT, double muT, double mB, double muB,
 double mC, double muC, double muLambda, double Q, double muH, double muJ, double muS,
 double R, double mu, int c[], long len, double lambda, double R0, double mu0, double delta0,
 double h, double tau, double tau2){
  double res;

f90singulardiffpiece_(hard, shape, gap, space, cum, &orderAlpha, &runAlpha, &order, &run, &nf,
&j3, &s3, &G3, &mZ, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &muLambda, &Q, &muH, &muJ,
&muS, &R, &mu, c, &lambda, &R0, &mu0, &delta0, &h, &tau, &tau2, &res);

return res;

}

extern double f90anomdim_(char const* str, int* nf, double* G4, double* result);

static void anomdim(char const* str, int nf, double G4){
  double result[4];

   f90anomdim_(str, &nf, &G4, result);

   MLPutRealList(stdlink, result, 4);

   MLEndPacket(stdlink);
}

extern double f90scoef_(char const* str, int* nf, double* result);

static void scoef(char const* str, int nf){
  double result[3];

   f90scoef_(str, &nf, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90fomass_(char const* shape, char const* current, double* m, double* Q,
          double* Mz, double* gammaZ, double* sin2ThetaW, double* tau, double* result);

static double fomass(char const* shape, char const* current, double m, double Q,
                     double Mz, double gammaZ, double sin2ThetaW, double tau){
  double result;

  f90fomass_(shape, current, &m, &Q, &Mz, &gammaZ, &sin2ThetaW, &tau, &result);
   return result;
}

extern double f90thrustns1loop_(double* tau, double* result);

static void thrustns1loop(double tau){
  double result[3];

   f90thrustns1loop_(&tau, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90thrustns2loop_(double *er, double* tau, double* result);

static void thrustns2loop(double er, double tau){
  double result[2];

   f90thrustns2loop_(&er, &tau, result);

   MLPutRealList(stdlink, result, 2);

   MLEndPacket(stdlink);
}

extern double f90polylog_(int* nf, double* z, double* result);

static double polylog(int n, double z){
  double res;

   f90polylog_(&n, &z, &res);

   return res;
}

extern double f90dilog_(double* z, double* result);

static double dilog(double z){
  double res;

   f90dilog_(&z, &res);

   return res;
}

extern double f90pfq_(double* a, int* clena, double* b, int* clenb, double* z, double* res);

static double pfq(double a[], long lena, double b[], long lenb, double z){
  int clena = lena; int clenb = lenb;
  double res;

   f90pfq_(a, &clena, b, &clenb, &z, &res);

   return res;
}

extern double f90elliptic3_(double* psi, double* k, double* c, double* result);

static double elliptic3(double psi, double k, double c){
  double res;

   f90elliptic3_(&psi, &k, &c, &res);

   return res;
}

extern double f90nglfunction_(int* nf, double* z, double* result);

static double nglfunction(int n, double z){
  double res;

   f90nglfunction_(&n, &z, &res);

   return res;
}

extern double f90mctop_(char const* str, double* mt, double* Q, int* n, int* k,
  double* x, double* result);

static double mctop(char const* str, double mt, double Q, int n, int k, double x){
  double res;

   f90mctop_(str, &mt, &Q, &n, &k, &x, &res);

   return res;
}


extern double f90breitunstable_(char const* str, double* mt, double* Q,
  double * gamma, int* n, int* k, double* x, double* result);

static double breitunstable(char const* str, double mt, double Q, double gamma,
  int n, int k, double x){
  double res;

   f90breitunstable_(str, &mt, &Q, &gamma, &n, &k, &x, &res);

   return res;
}

extern double f90deltamctop_(char const* str, double* mt, double* Q, double* result);

static double deltamctop(char const* str, double mt, double Q){
  double res;

   f90deltamctop_(str, &mt, &Q, &res);

   return res;
}

extern double f90complexpolylog_(int* nf, double* z1, double* z2, double* result);

static void complexpolylog(int n, double z1, double z2){
  double res[2];

   f90complexpolylog_(&n, &z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90nglsoft_(int* nf, double* z1, double* z2, double* result);

static void nglsoft(int n, double z1, double z2){
  double res[2];

   f90nglsoft_(&n, &z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90cli2_(double* z1, double* z2, double* result);

static void cli2(double z1, double z2){
  double res[2];

   f90cli2_(&z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90cli3_(double* z1, double* z2, double* result);

static void cli3(double z1, double z2){
  double res[2];

   f90cli3_(&z1, &z2, res);

   MLPutFunction(stdlink, "Complex", 2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);

   MLEndPacket(stdlink);
}

extern double f90deltagap_(char const* str, int* orderAlpha, int* runAlpha, int* runMass,
int* nf, double* mZ, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double* R, double* result);

static void deltagap(char const* str, int orderAlpha, int runAlpha, int runMass, int nf,
double mZ, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double mu, double R){
  double result[3];

   f90deltagap_(str, &orderAlpha, &runAlpha, &runMass, &nf, &mZ, &aMz, &mT, &muT, &mB,
   &muB, &mC, &muC, &mu, &R, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90delta_(char const* str, int* nf, double* mu, double* R, double* result);

static void delta(char const* str, int nf, double mu, double R){
  double result[3];

   f90delta_(str, &nf, &mu, &R, result);

   MLPutRealList(stdlink, result, 3);

   MLEndPacket(stdlink);
}

extern double f90diffdeltagap_(char const* str, char const* scheme, int* order, double*R0,
double* R1, double* mu0, double* mu1, double* muLambda, int* orderAlpha, int* runAlpha,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* result);

static double diffdeltagap(char const* str, char const* scheme, int order, double R0,
double R1, double mu0, double mu1, double muLambda, int orderAlpha, int runAlpha, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC){

  double result;

  f90diffdeltagap_(str, scheme, &order, &R0, &R1, &mu0, &mu1, &muLambda, &orderAlpha,
   &runAlpha, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &result);

   return result;
}

extern double f90diffdeltagapmass_(char const* str, int* order, double* R0, double* R1,
double* mu0, double* mu1, double* muM, double* muLambda1, double* muLambda2,
int* orderAlpha, int* runAlpha, int* runMass, int* nf, double* Mz, double* aMz, double* mT,
double* muT, double* mB, double* muB, double* mC, double* muC, double* result);

static double diffdeltagapmass(char const* str, int order, double R0, double R1,
double mu0, double mu1, double muM, double muLambda1, double muLambda2, int orderAlpha,
int runAlpha, int runMass, int nf, double Mz, double aMz, double mT, double muT, double mB,
double muB, double mC, double muC){

  double result;

  f90diffdeltagapmass_(str, &order, &R0, &R1, &mu0, &mu1, &muM, &muLambda1, &muLambda2,
  &orderAlpha, &runAlpha, &runMass, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
  &result);

   return result;
}

extern double f90coefmat_(char const* str, int* nf, double* s3, double* result);

static void coefmat(char const* str, int nf, double s3){
  double result[15];

   f90coefmat_(str, &nf, &s3, result);

   MLPutFunction(stdlink, "Partition", 2);
   MLPutRealList(stdlink, result, 15);
   MLPutInteger(stdlink, 5);

   MLEndPacket(stdlink);
}

extern double f90model_(double* c, int* clen, double* lambda, int* k, double* l,
  double* result);

static double model(double c[], long clen, double lambda, int k, double l){
  double res;
  int len = clen;

   f90model_(c, &len, &lambda, &k, &l, &res);

   return res;
}

extern double f90modelunstable_(char const* shape, double* mt, double* Q,
double* c, int* clen, double* lambda, int* n, int* k, double* l, double* result);

static double modelunstable(char const* shape, double mt, double Q, double c[], long clen,
double lambda, int n, int k, double l){
  double res;  int len = clen;

   f90modelunstable_(shape, &mt, &Q, c, &len, &lambda, &n, &k, &l, &res);

   return res;
}

extern double f90breitmodelunstable_(char const* shape, double* mt, double* Q,
double *gamma, double* c, int* clen, double* lambda, int* n, int* k, double* l,
double* result);

static double breitmodelunstable(char const* shape, double mt, double Q, double gamma,
  double c[], long clen, double lambda, int n, int k, double l){
  double res;  int len = clen;

   f90breitmodelunstable_(shape, &mt, &Q, &gamma, c, &len, &lambda, &n, &k, &l, &res);

   return res;
}

extern double f90modelunstablediff_(char const* shape, double* mt, double* Q,
double* c, int* clen, double* lambda, int* n, int* k, double* l, double* l2, double* result);

static double modelunstablediff(char const* shape, double mt, double Q, double c[], long clen,
double lambda, int n, int k, double l, double l2){
  double res;  int len = clen;

   f90modelunstablediff_(shape, &mt, &Q, c, &len, &lambda, &n, &k, &l, &l2, &res);

   return res;
}

extern double f90modeldiff_(double* c, int* clen, double* lambda, int* k, double* l,
double* l2, double* result);

static double modeldiff(double c[], long clen, double lambda, int k, double l, double l2){
  double res;
  int len = clen;

   f90modeldiff_(c, &len, &lambda, &k, &l, &l2, &res);

   return res;
}

extern double f90breitmodel_(double* c, int* clen, double* lambda, double* width,
int* k, double* l, double* result);

static double breitmodel(double c[], long clen, double lambda, double width, int k, double l){
  double res;
  int len = clen;

   f90breitmodel_(c, &len, &lambda, &width, &k, &l, &res);

   return res;
}

extern double f90breitmodeldiff_(double* c, int* clen, double* lambda, double* width,
int* k, double* l, double* l2, double* result);

static double breitmodeldiff(double c[], long clen, double lambda, double width,
int k, double l, double l2){
  double res;
  int len = clen;

   f90breitmodeldiff_(c, &len, &lambda, &width, &k, &l, &l2, &res);

   return res;
}

extern double f90taylor_(double* c, int* clen, double* lambda, int* k, double* result);

static double taylor(double c[], long clen, double lambda, int k){
  double res;
  int len = clen;

   f90taylor_(c, &len, &lambda, &k, &res);

   return res;
}

extern double f90momentmodel_(double* c, int* clen, double* lambda, int* k,
                              double* result);

static double momentmodel(double c[], long clen, double lambda, int k){
  double res;
  int len = clen;

   f90momentmodel_(c, &len, &lambda, &k, &res);

   return res;
}

extern double f90modelpiece_(int* c, double* lambda, int* k, double* l, double* result);

static double modelpiece(int c[], long clen, double lambda, int k, double l){
  double res;

   f90modelpiece_(c, &lambda, &k, &l, &res);

   return res;
}

extern double f90taylorpiece_(int* c, double* lambda, int* k, double* result);

static double taylorpiece(int c[], long clen, double lambda, int k){
  double res;

   f90taylorpiece_(c, &lambda, &k, &res);

   return res;
}

extern double f90wtilde_(int* order, int* nf, double* gamma, double* a0, double* a1,
                         double* res);

static double wtilde(int order, int nf, double gamma[], long ngamma, double a0, double a1){
  double res;

   f90wtilde_(&order, &nf, gamma, &a0, &a1, &res);

   return res;
}

extern double f90ktilde_(int* order, int* nf, double* gamma, double* a0, double* a1,
                         double* res);

static double ktilde(int order, int nf, double gamma[], long ngamma, double a0,
double a1){
  double res;

   f90ktilde_(&order, &nf, gamma, &a0, &a1, &res);

   return res;
}

extern double f90alphaqcd_(char const* str, char const* method, int* order, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static double alphaqcd(char const* str, char const* method, int order, int run, int nf, double Mz, double
aMz, double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res;

   f90alphaqcd_(str, method, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
   &mu, &res);

  return res;
}

extern double f90alphacomplex_(char const* str, char const* method, int* order,
int* run, int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB,
double* muB, double* mC, double* muC, double* muR, double* muI, double* res);

static void alphacomplex(char const* str, char const* method, int order, int run,
int nf, double Mz, double aMz, double mT, double muT, double mB, double muB,
double mC, double muC, double muR, double muI){

   double res[2];

   f90alphacomplex_(str, method, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
   &muB, &mC, &muC, &muR, &muI, res);

  //  MLPutRealList(stdlink, res, 2);
   MLPutFunction(stdlink,"Complex",2);
   MLPutReal(stdlink, res[0]);
   MLPutReal(stdlink, res[1]);
   MLEndPacket(stdlink);
}

extern double f90lambdaqcd_(char const* str, int* order, int* runAlpha, int* run, int* nf,
double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mu, double* res);

static double lambdaqcd(char const* str, int order, int runAlpha, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double mu){

   double res;

   f90lambdaqcd_(str, &order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC,
   &muC, &mu, &res);

  return res;
}

extern double f90msbarmass_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static double msbarmass(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res;

   f90msbarmass_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC,
                 &muC, &mu, &res);

  return res;
}

extern double f90polemass_(int* orderAlpha, int* runAlpha, int*order, int* run, int* nf,
double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mu, double* res);

static double polemass(int orderAlpha, int runAlpha, int order, int run, int nf, double Mz,
double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double mu){

   double res;

   f90polemass_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
                 &muB, &mC, &muC, &mu, &res);

  return res;
}

extern double f90msbarmasslow_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static double msbarmasslow(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res;

   f90msbarmasslow_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC,
                 &muC, &mu, &res);

  return res;
}

extern double f90msrmass_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* R, double* res);

static double msrmass(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu, double R){

   double res;

   f90msrmass_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
               &mu, &R, &res);

  return res;
}

extern double f90mmfrommsr_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* R, double* res);

static double mmfrommsr(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu, double R){

   double res;

   f90mmfrommsr_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC, &muC,
               &mu, &R, &res);

  return res;
}

extern double f90mmfrommsrnatural_(int* orderAlpha, int* runAlpha, int* order, int* run,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double* R, double* res);

static double mmfrommsrnatural(int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double mu, double R){

   double res;

   f90mmfrommsrnatural_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT,
   &mB, &muB, &mC, &muC, &mu, &R, &res);

  return res;
}

extern double f90msrnaturalmass_(int* orderAlpha, int* runAlpha, int* order, int* run,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double* R, double* res);

static double msrnaturalmass(int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double mu, double R){

   double res;

   f90msrnaturalmass_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT,
   &mB, &muB, &mC, &muC, &mu, &R, &res);

  return res;
}

extern double f90jetmass_(int* orderAlpha, int* runAlpha, int* order, int* run, int* nf,
double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* muLambda, double* R, double* mu, double* res);

static double jetmass(int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double muLambda, double R, double mu){

   double res;

   f90jetmass_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
               &muB, &mC, &muC, &muLambda, &R, &mu, &res);

  return res;
}

extern double f90mmfromjetmass_(int* orderAlpha, int* runAlpha, int* order, int* run,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* muLambda, double* R, double* mu, double* res);

static double mmfromjetmass(int orderAlpha, int runAlpha, int order, int run, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC,
double muC, double muLambda, double R, double mu){

   double res;

   f90mmfromjetmass_(&orderAlpha, &runAlpha, &order, &run, &nf, &Mz, &aMz, &mT, &muT, &mB,
               &muB, &mC, &muC, &muLambda, &R, &mu, &res);

  return res;
}

extern double f90deltamsbar_(int* order, int* runAlpha, int* run, int* nf, double* Mz,
double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC, double* muC,
double* mu, double* res);

static void deltamsbar(int order, int runAlpha, int run, int nf, double Mz, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu){

   double res[3];

   f90deltamsbar_(&order, &runAlpha, &run, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB, &mC,
   &muC, &mu, res);

   MLPutRealList(stdlink, res, 3);

   MLEndPacket(stdlink);
}

extern double f90rhad_(char const* str, int* orderAlpha, int* runAlpha, int* order,
int* nf, double* Mz, double* aMz, double* mT, double* muT, double* mB, double* muB,
double* mC, double* muC, double* mu, double * Q, double* res);

static double rhad(char const* str, int orderAlpha, int runAlpha, int order, int nf,
double Mz, double aMz, double mT, double muT, double mB, double muB, double mC, double muC,
double mu, double Q){

   double res;

   f90rhad_(str, &orderAlpha, &runAlpha, &order, &nf, &Mz, &aMz, &mT, &muT, &mB, &muB,
   &mC, &muC, &mu, &Q, &res);

  return res;
}

extern double f90rhadmass_(char const* str, char const* curr, int* orderAlpha,
int* runAlpha, int* runMass, int* order, int* nf, double* Mz, double* gammaZ,
double* sinW, double* aMz, double* mT, double* muT, double* mB, double* muB, double* mC,
double* muC, double* mu, double * Q, double* res);

static double rhadmass(char const* str, char const* curr, int orderAlpha, int runAlpha,
int runMass, int order, int nf, double Mz, double gammaZ, double sinW, double aMz,
double mT, double muT, double mB, double muB, double mC, double muC, double mu, double Q){

   double res;

   f90rhadmass_(str, curr, &orderAlpha, &runAlpha, &runMass, &order, &nf, &Mz, &gammaZ,
                &sinW, &aMz, &mT, &muT, &mB, &muB, &mC, &muC, &mu, &Q, &res);

  return res;
}

extern double f90hyper2f1_(double *a, double *b, double *c, double *x, double *res);

static double hyper2f1( double a, double b, double c, double x){
  double res;

   f90hyper2f1_(&a, &b, &c, &x, &res);

   return res;
}

extern double f90intecorre_(double *b, double *x0, double *x1, double *res);

static double intecorre( double b, double x0, double x1){
  double res;

   f90intecorre_(&b, &x0, &x1, &res);

   return res;
}

extern double f90hyperf32exact_(double *w, double *x, double *res);

static double hyperf32exact(double w, double x){
  double res;

   f90hyperf32exact_(&w, &x, &res);

   return res;
}

int main(int argc, char *argv[]){
    return MLMain(argc, argv);
}
# line 2558 "/Users/vicent/GitHubProjects/Caliper/src/Caliper.tm.c"


void ewfactors P(( int _tp1, double _tp2, double _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr0( MLINK mlp)
#else
static int _tr0(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	ewfactors(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr0 */


void legendrelist P(( int _tp1, int _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr1( MLINK mlp)
#else
static int _tr1(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double _tp3;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	legendrelist(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr1 */


void qlegendrelist P(( int _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr2( MLINK mlp)
#else
static int _tr2(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	qlegendrelist(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr2 */


void gammar P(( const char * _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr3( MLINK mlp)
#else
static int _tr3(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	gammar(_tp1, _tp2);

	res = 1;
L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr3 */


double masslessprof P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, int _tp39, double * _tp40, long _tpl40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46));

#if MLPROTOTYPES
static int _tr4( MLINK mlp)
#else
static int _tr4(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	int _tp39;
	double * _tp40;
	long _tpl40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetRealList( mlp, &_tp40, &_tpl40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLNewPacket(mlp) ) goto L46;

	_rp0 = masslessprof(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tpl40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L46: L45: L44: L43: L42: L41: L40:	MLReleaseReal64List(mlp, _tp40, _tpl40);
L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr4 */


double findorigin P(( const char * _tp1, const char * _tp2, int _tp3, int _tp4, int _tp5, int _tp6, int _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33));

#if MLPROTOTYPES
static int _tr5( MLINK mlp)
#else
static int _tr5(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	int _tp6;
	int _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLNewPacket(mlp) ) goto L33;

	_rp0 = findorigin(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr5 */


void masslessprofpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, int _tp38, int _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45));

#if MLPROTOTYPES
static int _tr6( MLINK mlp)
#else
static int _tr6(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	int _tp38;
	int _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetInteger( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLNewPacket(mlp) ) goto L45;

	masslessprofpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45);

	res = 1;
L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr6 */


void masslessprofpiecelist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, int _tp38, int _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double * _tp45, long _tpl45));

#if MLPROTOTYPES
static int _tr7( MLINK mlp)
#else
static int _tr7(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	int _tp38;
	int _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double * _tp45;
	long _tpl45;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetInteger( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetRealList( mlp, &_tp45, &_tpl45) ) goto L44;
	if ( ! MLNewPacket(mlp) ) goto L45;

	masslessprofpiecelist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tpl45);

	res = 1;
L45:	MLReleaseReal64List(mlp, _tp45, _tpl45);
L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr7 */


void masslesspiecebin P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, int _tp38, int _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double * _tp45, long _tpl45, int _tp46));

#if MLPROTOTYPES
static int _tr8( MLINK mlp)
#else
static int _tr8(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	int _tp38;
	int _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double * _tp45;
	long _tpl45;
	int _tp46;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetInteger( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetRealList( mlp, &_tp45, &_tpl45) ) goto L44;
	if ( ! MLGetInteger( mlp, &_tp46) ) goto L45;
	if ( ! MLNewPacket(mlp) ) goto L46;

	masslesspiecebin(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tpl45, _tp46);

	res = 1;
L46: L45:	MLReleaseReal64List(mlp, _tp45, _tpl45);
L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr8 */


void masslessprofdiffpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, int _tp38, int _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46));

#if MLPROTOTYPES
static int _tr9( MLINK mlp)
#else
static int _tr9(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	int _tp38;
	int _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetInteger( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLNewPacket(mlp) ) goto L46;

	masslessprofdiffpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46);

	res = 1;
L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr9 */


void masslessproflist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, int _tp39, double * _tp40, long _tpl40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double * _tp46, long _tpl46));

#if MLPROTOTYPES
static int _tr10( MLINK mlp)
#else
static int _tr10(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	int _tp39;
	double * _tp40;
	long _tpl40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double * _tp46;
	long _tpl46;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetRealList( mlp, &_tp40, &_tpl40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetRealList( mlp, &_tp46, &_tpl46) ) goto L45;
	if ( ! MLNewPacket(mlp) ) goto L46;

	masslessproflist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tpl40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tpl46);

	res = 1;
L46:	MLReleaseReal64List(mlp, _tp46, _tpl46);
L45: L44: L43: L42: L41: L40:	MLReleaseReal64List(mlp, _tp40, _tpl40);
L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr10 */


void masslessbinlist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, int _tp39, double * _tp40, long _tpl40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double * _tp46, long _tpl46, int _tp47));

#if MLPROTOTYPES
static int _tr11( MLINK mlp)
#else
static int _tr11(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	int _tp39;
	double * _tp40;
	long _tpl40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double * _tp46;
	long _tpl46;
	int _tp47;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetRealList( mlp, &_tp40, &_tpl40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetRealList( mlp, &_tp46, &_tpl46) ) goto L45;
	if ( ! MLGetInteger( mlp, &_tp47) ) goto L46;
	if ( ! MLNewPacket(mlp) ) goto L47;

	masslessbinlist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tpl40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tpl46, _tp47);

	res = 1;
L47: L46:	MLReleaseReal64List(mlp, _tp46, _tpl46);
L45: L44: L43: L42: L41: L40:	MLReleaseReal64List(mlp, _tp40, _tpl40);
L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr11 */


double masslessprofdiff P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, int _tp39, double * _tp40, long _tpl40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47));

#if MLPROTOTYPES
static int _tr12( MLINK mlp)
#else
static int _tr12(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	int _tp39;
	double * _tp40;
	long _tpl40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetInteger( mlp, &_tp39) ) goto L38;
	if ( ! MLGetRealList( mlp, &_tp40, &_tpl40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLNewPacket(mlp) ) goto L47;

	_rp0 = masslessprofdiff(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tpl40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L47: L46: L45: L44: L43: L42: L41: L40:	MLReleaseReal64List(mlp, _tp40, _tpl40);
L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr12 */


double masslessmoment P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, int _tp38, double * _tp39, long _tpl39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, int _tp47));

#if MLPROTOTYPES
static int _tr13( MLINK mlp)
#else
static int _tr13(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	int _tp38;
	double * _tp39;
	long _tpl39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	int _tp47;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetInteger( mlp, &_tp38) ) goto L37;
	if ( ! MLGetRealList( mlp, &_tp39, &_tpl39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetInteger( mlp, &_tp47) ) goto L46;
	if ( ! MLNewPacket(mlp) ) goto L47;

	_rp0 = masslessmoment(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tpl39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L47: L46: L45: L44: L43: L42: L41: L40: L39:	MLReleaseReal64List(mlp, _tp39, _tpl39);
L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr13 */


void profiles P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, int _tp15, double _tp16));

#if MLPROTOTYPES
static int _tr14( MLINK mlp)
#else
static int _tr14(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	int _tp15;
	double _tp16;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLNewPacket(mlp) ) goto L16;

	profiles(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16);

	res = 1;
L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr14 */


void profilesmass P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, int _tp19, const char * _tp20, const char * _tp21, double _tp22));

#if MLPROTOTYPES
static int _tr15( MLINK mlp)
#else
static int _tr15(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	int _tp19;
	const char * _tp20;
	const char * _tp21;
	double _tp22;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetString( mlp, &_tp20) ) goto L19;
	if ( ! MLGetString( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLNewPacket(mlp) ) goto L22;

	profilesmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22);

	res = 1;
L22: L21:	MLReleaseString(mlp, _tp21);
L20:	MLReleaseString(mlp, _tp20);
L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr15 */


double mctop P(( const char * _tp1, double _tp2, double _tp3, int _tp4, int _tp5, double _tp6));

#if MLPROTOTYPES
static int _tr16( MLINK mlp)
#else
static int _tr16(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	double _tp2;
	double _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLNewPacket(mlp) ) goto L6;

	_rp0 = mctop(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L6: L5: L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr16 */


double breitunstable P(( const char * _tp1, double _tp2, double _tp3, double _tp4, int _tp5, int _tp6, double _tp7));

#if MLPROTOTYPES
static int _tr17( MLINK mlp)
#else
static int _tr17(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	int _tp5;
	int _tp6;
	double _tp7;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLNewPacket(mlp) ) goto L7;

	_rp0 = breitunstable(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L7: L6: L5: L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr17 */


double deltamctop P(( const char * _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr18( MLINK mlp)
#else
static int _tr18(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	double _tp2;
	double _tp3;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	_rp0 = deltamctop(_tp1, _tp2, _tp3);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr18 */


void massintersection P(( double _tp1, double _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, const char * _tp18, const char * _tp19));

#if MLPROTOTYPES
static int _tr19( MLINK mlp)
#else
static int _tr19(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	const char * _tp18;
	const char * _tp19;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetString( mlp, &_tp18) ) goto L17;
	if ( ! MLGetString( mlp, &_tp19) ) goto L18;
	if ( ! MLNewPacket(mlp) ) goto L19;

	massintersection(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19);

	res = 1;
L19:	MLReleaseString(mlp, _tp19);
L18:	MLReleaseString(mlp, _tp18);
L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr19 */


double nsmasspiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, int * _tp32, long _tpl32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40));

#if MLPROTOTYPES
static int _tr20( MLINK mlp)
#else
static int _tr20(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	int * _tp32;
	long _tpl32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetIntegerList( mlp, &_tp32, &_tpl32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLNewPacket(mlp) ) goto L40;

	_rp0 = nsmasspiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tpl32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L40: L39: L38: L37: L36: L35: L34: L33: L32:	MLReleaseInteger32List( mlp, _tp32, _tpl32);
L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr20 */


double nsmassdiffpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, int * _tp32, long _tpl32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41));

#if MLPROTOTYPES
static int _tr21( MLINK mlp)
#else
static int _tr21(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	int * _tp32;
	long _tpl32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetIntegerList( mlp, &_tp32, &_tpl32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLNewPacket(mlp) ) goto L41;

	_rp0 = nsmassdiffpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tpl32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L41: L40: L39: L38: L37: L36: L35: L34: L33: L32:	MLReleaseInteger32List( mlp, _tp32, _tpl32);
L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr21 */


double nsmass P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, int _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double * _tp33, long _tpl33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41));

#if MLPROTOTYPES
static int _tr22( MLINK mlp)
#else
static int _tr22(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	int _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double * _tp33;
	long _tpl33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetRealList( mlp, &_tp33, &_tpl33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLNewPacket(mlp) ) goto L41;

	_rp0 = nsmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tpl33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L41: L40: L39: L38: L37: L36: L35: L34: L33:	MLReleaseReal64List(mlp, _tp33, _tpl33);
L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr22 */


double nsmassdiff P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, int _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double * _tp33, long _tpl33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42));

#if MLPROTOTYPES
static int _tr23( MLINK mlp)
#else
static int _tr23(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	int _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double * _tp33;
	long _tpl33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetRealList( mlp, &_tp33, &_tpl33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLNewPacket(mlp) ) goto L42;

	_rp0 = nsmassdiff(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tpl33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L42: L41: L40: L39: L38: L37: L36: L35: L34: L33:	MLReleaseReal64List(mlp, _tp33, _tpl33);
L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr23 */


double hjmnsmass P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double * _tp32, long _tpl32, int _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41));

#if MLPROTOTYPES
static int _tr24( MLINK mlp)
#else
static int _tr24(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double * _tp32;
	long _tpl32;
	int _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetRealList( mlp, &_tp32, &_tpl32) ) goto L31;
	if ( ! MLGetInteger( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLNewPacket(mlp) ) goto L41;

	_rp0 = hjmnsmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tpl32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L41: L40: L39: L38: L37: L36: L35: L34: L33: L32:	MLReleaseReal64List(mlp, _tp32, _tpl32);
L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr24 */


double hjmnsmassdiff P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double * _tp32, long _tpl32, int _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42));

#if MLPROTOTYPES
static int _tr25( MLINK mlp)
#else
static int _tr25(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double * _tp32;
	long _tpl32;
	int _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetRealList( mlp, &_tp32, &_tpl32) ) goto L31;
	if ( ! MLGetInteger( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLNewPacket(mlp) ) goto L42;

	_rp0 = hjmnsmassdiff(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tpl32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32:	MLReleaseReal64List(mlp, _tp32, _tpl32);
L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr25 */


double singular P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double * _tp30, long _tpl30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36));

#if MLPROTOTYPES
static int _tr26( MLINK mlp)
#else
static int _tr26(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double * _tp30;
	long _tpl30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetRealList( mlp, &_tp30, &_tpl30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLNewPacket(mlp) ) goto L36;

	_rp0 = singular(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tpl30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L36: L35: L34: L33: L32: L31: L30:	MLReleaseReal64List(mlp, _tp30, _tpl30);
L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr26 */


void singularlist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, int _tp6, int _tp7, int _tp8, int _tp9, int _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, int _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35));

#if MLPROTOTYPES
static int _tr27( MLINK mlp)
#else
static int _tr27(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	int _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetInteger( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLNewPacket(mlp) ) goto L35;

	singularlist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35);

	res = 1;
L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr27 */


double singulardiff P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double * _tp30, long _tpl30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37));

#if MLPROTOTYPES
static int _tr28( MLINK mlp)
#else
static int _tr28(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double * _tp30;
	long _tpl30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetRealList( mlp, &_tp30, &_tpl30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLNewPacket(mlp) ) goto L37;

	_rp0 = singulardiff(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tpl30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L37: L36: L35: L34: L33: L32: L31: L30:	MLReleaseReal64List(mlp, _tp30, _tpl30);
L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr28 */


double massiveprof P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, double * _tp54, long _tpl54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double _tp62));

#if MLPROTOTYPES
static int _tr29( MLINK mlp)
#else
static int _tr29(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	double * _tp54;
	long _tpl54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double _tp62;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetRealList( mlp, &_tp54, &_tpl54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetReal( mlp, &_tp62) ) goto L61;
	if ( ! MLNewPacket(mlp) ) goto L62;

	_rp0 = massiveprof(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tpl54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L62: L61: L60: L59: L58: L57: L56: L55: L54:	MLReleaseReal64List(mlp, _tp54, _tpl54);
L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr29 */


double massorigin P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, int _tp5, int _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43));

#if MLPROTOTYPES
static int _tr30( MLINK mlp)
#else
static int _tr30(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	int _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLNewPacket(mlp) ) goto L43;

	_rp0 = massorigin(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr30 */


void massiveprofpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, int _tp54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double _tp62));

#if MLPROTOTYPES
static int _tr31( MLINK mlp)
#else
static int _tr31(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	int _tp54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double _tp62;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetInteger( mlp, &_tp54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetReal( mlp, &_tp62) ) goto L61;
	if ( ! MLNewPacket(mlp) ) goto L62;

	massiveprofpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62);

	res = 1;
L62: L61: L60: L59: L58: L57: L56: L55: L54: L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr31 */


void massiveprofpiecelist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, int _tp54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double * _tp62, long _tpl62));

#if MLPROTOTYPES
static int _tr32( MLINK mlp)
#else
static int _tr32(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	int _tp54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double * _tp62;
	long _tpl62;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetInteger( mlp, &_tp54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetRealList( mlp, &_tp62, &_tpl62) ) goto L61;
	if ( ! MLNewPacket(mlp) ) goto L62;

	massiveprofpiecelist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62, _tpl62);

	res = 1;
L62:	MLReleaseReal64List(mlp, _tp62, _tpl62);
L61: L60: L59: L58: L57: L56: L55: L54: L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr32 */


void massivepiecebin P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, int _tp54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double * _tp62, long _tpl62, int _tp63));

#if MLPROTOTYPES
static int _tr33( MLINK mlp)
#else
static int _tr33(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	int _tp54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double * _tp62;
	long _tpl62;
	int _tp63;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetInteger( mlp, &_tp54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetRealList( mlp, &_tp62, &_tpl62) ) goto L61;
	if ( ! MLGetInteger( mlp, &_tp63) ) goto L62;
	if ( ! MLNewPacket(mlp) ) goto L63;

	massivepiecebin(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62, _tpl62, _tp63);

	res = 1;
L63: L62:	MLReleaseReal64List(mlp, _tp62, _tpl62);
L61: L60: L59: L58: L57: L56: L55: L54: L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr33 */


void massiveprofdiffpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, int _tp54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double _tp62, double _tp63));

#if MLPROTOTYPES
static int _tr34( MLINK mlp)
#else
static int _tr34(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	int _tp54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double _tp62;
	double _tp63;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetInteger( mlp, &_tp54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetReal( mlp, &_tp62) ) goto L61;
	if ( ! MLGetReal( mlp, &_tp63) ) goto L62;
	if ( ! MLNewPacket(mlp) ) goto L63;

	massiveprofdiffpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62, _tp63);

	res = 1;
L63: L62: L61: L60: L59: L58: L57: L56: L55: L54: L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr34 */


double massiveprofdiff P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, double * _tp54, long _tpl54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double _tp62, double _tp63));

#if MLPROTOTYPES
static int _tr35( MLINK mlp)
#else
static int _tr35(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	double * _tp54;
	long _tpl54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double _tp62;
	double _tp63;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetRealList( mlp, &_tp54, &_tpl54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetReal( mlp, &_tp62) ) goto L61;
	if ( ! MLGetReal( mlp, &_tp63) ) goto L62;
	if ( ! MLNewPacket(mlp) ) goto L63;

	_rp0 = massiveprofdiff(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tpl54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62, _tp63);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L63: L62: L61: L60: L59: L58: L57: L56: L55: L54:	MLReleaseReal64List(mlp, _tp54, _tpl54);
L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr35 */


void massiveproflist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, double * _tp54, long _tpl54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double * _tp62, long _tpl62));

#if MLPROTOTYPES
static int _tr36( MLINK mlp)
#else
static int _tr36(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	double * _tp54;
	long _tpl54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double * _tp62;
	long _tpl62;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetRealList( mlp, &_tp54, &_tpl54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetRealList( mlp, &_tp62, &_tpl62) ) goto L61;
	if ( ! MLNewPacket(mlp) ) goto L62;

	massiveproflist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tpl54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62, _tpl62);

	res = 1;
L62:	MLReleaseReal64List(mlp, _tp62, _tpl62);
L61: L60: L59: L58: L57: L56: L55: L54:	MLReleaseReal64List(mlp, _tp54, _tpl54);
L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr36 */


void massivebinlist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, const char * _tp11, double _tp12, double _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, int _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51, int _tp52, double _tp53, double * _tp54, long _tpl54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double * _tp62, long _tpl62, int _tp63));

#if MLPROTOTYPES
static int _tr37( MLINK mlp)
#else
static int _tr37(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	const char * _tp11;
	double _tp12;
	double _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	int _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	int _tp52;
	double _tp53;
	double * _tp54;
	long _tpl54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double * _tp62;
	long _tpl62;
	int _tp63;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetString( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetInteger( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLGetInteger( mlp, &_tp52) ) goto L51;
	if ( ! MLGetReal( mlp, &_tp53) ) goto L52;
	if ( ! MLGetRealList( mlp, &_tp54, &_tpl54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetRealList( mlp, &_tp62, &_tpl62) ) goto L61;
	if ( ! MLGetInteger( mlp, &_tp63) ) goto L62;
	if ( ! MLNewPacket(mlp) ) goto L63;

	massivebinlist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tp54, _tpl54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62, _tpl62, _tp63);

	res = 1;
L63: L62:	MLReleaseReal64List(mlp, _tp62, _tpl62);
L61: L60: L59: L58: L57: L56: L55: L54:	MLReleaseReal64List(mlp, _tp54, _tpl54);
L53: L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11:	MLReleaseString(mlp, _tp11);
L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr37 */


double massivemoment P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, double _tp11, double _tp12, int _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, int _tp51, double _tp52, double * _tp53, long _tpl53, double _tp54, double _tp55, double _tp56, double _tp57, double _tp58, double _tp59, double _tp60, double _tp61, double _tp62, int _tp63));

#if MLPROTOTYPES
static int _tr38( MLINK mlp)
#else
static int _tr38(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	double _tp11;
	double _tp12;
	int _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	int _tp51;
	double _tp52;
	double * _tp53;
	long _tpl53;
	double _tp54;
	double _tp55;
	double _tp56;
	double _tp57;
	double _tp58;
	double _tp59;
	double _tp60;
	double _tp61;
	double _tp62;
	int _tp63;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetInteger( mlp, &_tp51) ) goto L50;
	if ( ! MLGetReal( mlp, &_tp52) ) goto L51;
	if ( ! MLGetRealList( mlp, &_tp53, &_tpl53) ) goto L52;
	if ( ! MLGetReal( mlp, &_tp54) ) goto L53;
	if ( ! MLGetReal( mlp, &_tp55) ) goto L54;
	if ( ! MLGetReal( mlp, &_tp56) ) goto L55;
	if ( ! MLGetReal( mlp, &_tp57) ) goto L56;
	if ( ! MLGetReal( mlp, &_tp58) ) goto L57;
	if ( ! MLGetReal( mlp, &_tp59) ) goto L58;
	if ( ! MLGetReal( mlp, &_tp60) ) goto L59;
	if ( ! MLGetReal( mlp, &_tp61) ) goto L60;
	if ( ! MLGetReal( mlp, &_tp62) ) goto L61;
	if ( ! MLGetInteger( mlp, &_tp63) ) goto L62;
	if ( ! MLNewPacket(mlp) ) goto L63;

	_rp0 = massivemoment(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51, _tp52, _tp53, _tpl53, _tp54, _tp55, _tp56, _tp57, _tp58, _tp59, _tp60, _tp61, _tp62, _tp63);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L63: L62: L61: L60: L59: L58: L57: L56: L55: L54: L53:	MLReleaseReal64List(mlp, _tp53, _tpl53);
L52: L51: L50: L49: L48: L47: L46: L45: L44: L43: L42: L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr38 */


double singularmass P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, double _tp11, double _tp12, int _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double * _tp42, long _tpl42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50));

#if MLPROTOTYPES
static int _tr39( MLINK mlp)
#else
static int _tr39(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	double _tp11;
	double _tp12;
	int _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double * _tp42;
	long _tpl42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetRealList( mlp, &_tp42, &_tpl42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLNewPacket(mlp) ) goto L50;

	_rp0 = singularmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tpl42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L50: L49: L48: L47: L46: L45: L44: L43: L42:	MLReleaseReal64List(mlp, _tp42, _tpl42);
L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr39 */


double singularmassdiff P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, double _tp11, double _tp12, int _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double * _tp42, long _tpl42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51));

#if MLPROTOTYPES
static int _tr40( MLINK mlp)
#else
static int _tr40(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	double _tp11;
	double _tp12;
	int _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double * _tp42;
	long _tpl42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetRealList( mlp, &_tp42, &_tpl42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLNewPacket(mlp) ) goto L51;

	_rp0 = singularmassdiff(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tpl42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L51: L50: L49: L48: L47: L46: L45: L44: L43: L42:	MLReleaseReal64List(mlp, _tp42, _tpl42);
L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr40 */


double massnondist P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, int _tp14, int _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double * _tp35, long _tpl35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41));

#if MLPROTOTYPES
static int _tr41( MLINK mlp)
#else
static int _tr41(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	int _tp14;
	int _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double * _tp35;
	long _tpl35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetRealList( mlp, &_tp35, &_tpl35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLNewPacket(mlp) ) goto L41;

	_rp0 = massnondist(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tpl35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L41: L40: L39: L38: L37: L36: L35:	MLReleaseReal64List(mlp, _tp35, _tpl35);
L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr41 */


double massnondistdiff P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, int _tp14, int _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double * _tp35, long _tpl35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, double _tp42));

#if MLPROTOTYPES
static int _tr42( MLINK mlp)
#else
static int _tr42(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	int _tp14;
	int _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double * _tp35;
	long _tpl35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _tp42;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetRealList( mlp, &_tp35, &_tpl35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetReal( mlp, &_tp42) ) goto L41;
	if ( ! MLNewPacket(mlp) ) goto L42;

	_rp0 = massnondistdiff(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tpl35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L42: L41: L40: L39: L38: L37: L36: L35:	MLReleaseReal64List(mlp, _tp35, _tpl35);
L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr42 */


double massnondistpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, int _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, int * _tp34, long _tpl34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40));

#if MLPROTOTYPES
static int _tr43( MLINK mlp)
#else
static int _tr43(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	int _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	int * _tp34;
	long _tpl34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetIntegerList( mlp, &_tp34, &_tpl34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLNewPacket(mlp) ) goto L40;

	_rp0 = massnondistpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tpl34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L40: L39: L38: L37: L36: L35: L34:	MLReleaseInteger32List( mlp, _tp34, _tpl34);
L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr43 */


double massnondistdiffpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, int _tp13, int _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, int * _tp34, long _tpl34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41));

#if MLPROTOTYPES
static int _tr44( MLINK mlp)
#else
static int _tr44(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	int _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	int * _tp34;
	long _tpl34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetIntegerList( mlp, &_tp34, &_tpl34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLNewPacket(mlp) ) goto L41;

	_rp0 = massnondistdiffpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tpl34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L41: L40: L39: L38: L37: L36: L35: L34:	MLReleaseInteger32List( mlp, _tp34, _tpl34);
L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr44 */


double singularmasspiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, double _tp11, double _tp12, int _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, int * _tp42, long _tpl42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50));

#if MLPROTOTYPES
static int _tr45( MLINK mlp)
#else
static int _tr45(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	double _tp11;
	double _tp12;
	int _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	int * _tp42;
	long _tpl42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetIntegerList( mlp, &_tp42, &_tpl42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLNewPacket(mlp) ) goto L50;

	_rp0 = singularmasspiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tpl42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L50: L49: L48: L47: L46: L45: L44: L43: L42:	MLReleaseInteger32List( mlp, _tp42, _tpl42);
L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr45 */


double singularmassdiffpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, const char * _tp7, const char * _tp8, const char * _tp9, const char * _tp10, double _tp11, double _tp12, int _tp13, int _tp14, int _tp15, int _tp16, int _tp17, int _tp18, int _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41, int * _tp42, long _tpl42, double _tp43, double _tp44, double _tp45, double _tp46, double _tp47, double _tp48, double _tp49, double _tp50, double _tp51));

#if MLPROTOTYPES
static int _tr46( MLINK mlp)
#else
static int _tr46(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	const char * _tp7;
	const char * _tp8;
	const char * _tp9;
	const char * _tp10;
	double _tp11;
	double _tp12;
	int _tp13;
	int _tp14;
	int _tp15;
	int _tp16;
	int _tp17;
	int _tp18;
	int _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	int * _tp42;
	long _tpl42;
	double _tp43;
	double _tp44;
	double _tp45;
	double _tp46;
	double _tp47;
	double _tp48;
	double _tp49;
	double _tp50;
	double _tp51;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetString( mlp, &_tp7) ) goto L6;
	if ( ! MLGetString( mlp, &_tp8) ) goto L7;
	if ( ! MLGetString( mlp, &_tp9) ) goto L8;
	if ( ! MLGetString( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetInteger( mlp, &_tp14) ) goto L13;
	if ( ! MLGetInteger( mlp, &_tp15) ) goto L14;
	if ( ! MLGetInteger( mlp, &_tp16) ) goto L15;
	if ( ! MLGetInteger( mlp, &_tp17) ) goto L16;
	if ( ! MLGetInteger( mlp, &_tp18) ) goto L17;
	if ( ! MLGetInteger( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLGetIntegerList( mlp, &_tp42, &_tpl42) ) goto L41;
	if ( ! MLGetReal( mlp, &_tp43) ) goto L42;
	if ( ! MLGetReal( mlp, &_tp44) ) goto L43;
	if ( ! MLGetReal( mlp, &_tp45) ) goto L44;
	if ( ! MLGetReal( mlp, &_tp46) ) goto L45;
	if ( ! MLGetReal( mlp, &_tp47) ) goto L46;
	if ( ! MLGetReal( mlp, &_tp48) ) goto L47;
	if ( ! MLGetReal( mlp, &_tp49) ) goto L48;
	if ( ! MLGetReal( mlp, &_tp50) ) goto L49;
	if ( ! MLGetReal( mlp, &_tp51) ) goto L50;
	if ( ! MLNewPacket(mlp) ) goto L51;

	_rp0 = singularmassdiffpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41, _tp42, _tpl42, _tp43, _tp44, _tp45, _tp46, _tp47, _tp48, _tp49, _tp50, _tp51);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L51: L50: L49: L48: L47: L46: L45: L44: L43: L42:	MLReleaseInteger32List( mlp, _tp42, _tpl42);
L41: L40: L39: L38: L37: L36: L35: L34: L33: L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10:	MLReleaseString(mlp, _tp10);
L9:	MLReleaseString(mlp, _tp9);
L8:	MLReleaseString(mlp, _tp8);
L7:	MLReleaseString(mlp, _tp7);
L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr46 */


double singularhjm P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, int _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double * _tp32, long _tpl32, int _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39));

#if MLPROTOTYPES
static int _tr47( MLINK mlp)
#else
static int _tr47(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double * _tp32;
	long _tpl32;
	int _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetRealList( mlp, &_tp32, &_tpl32) ) goto L31;
	if ( ! MLGetInteger( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLNewPacket(mlp) ) goto L39;

	_rp0 = singularhjm(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tpl32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L39: L38: L37: L36: L35: L34: L33: L32:	MLReleaseReal64List(mlp, _tp32, _tpl32);
L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr47 */


double singularhjm1d P(( const char * _tp1, const char * _tp2, const char * _tp3, int _tp4, int _tp5, int _tp6, int _tp7, int _tp8, int _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double * _tp30, long _tpl30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36));

#if MLPROTOTYPES
static int _tr48( MLINK mlp)
#else
static int _tr48(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	int _tp4;
	int _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double * _tp30;
	long _tpl30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetRealList( mlp, &_tp30, &_tpl30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLNewPacket(mlp) ) goto L36;

	_rp0 = singularhjm1d(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tpl30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L36: L35: L34: L33: L32: L31: L30:	MLReleaseReal64List(mlp, _tp30, _tpl30);
L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr48 */


double singularhjm1dpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, int _tp4, int _tp5, int _tp6, int _tp7, int _tp8, int _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, int * _tp30, long _tpl30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36));

#if MLPROTOTYPES
static int _tr49( MLINK mlp)
#else
static int _tr49(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	int _tp4;
	int _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	int * _tp30;
	long _tpl30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetIntegerList( mlp, &_tp30, &_tpl30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLNewPacket(mlp) ) goto L36;

	_rp0 = singularhjm1dpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tpl30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L36: L35: L34: L33: L32: L31: L30:	MLReleaseInteger32List( mlp, _tp30, _tpl30);
L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr49 */


double singulardouble P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, const char * _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, int _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, double _tp32, double * _tp33, long _tpl33, int _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39, double _tp40, double _tp41));

#if MLPROTOTYPES
static int _tr50( MLINK mlp)
#else
static int _tr50(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	const char * _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	double _tp32;
	double * _tp33;
	long _tpl33;
	int _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _tp40;
	double _tp41;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetString( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetRealList( mlp, &_tp33, &_tpl33) ) goto L32;
	if ( ! MLGetInteger( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLGetReal( mlp, &_tp40) ) goto L39;
	if ( ! MLGetReal( mlp, &_tp41) ) goto L40;
	if ( ! MLNewPacket(mlp) ) goto L41;

	_rp0 = singulardouble(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tp33, _tpl33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39, _tp40, _tp41);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L41: L40: L39: L38: L37: L36: L35: L34: L33:	MLReleaseReal64List(mlp, _tp33, _tpl33);
L32: L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6:	MLReleaseString(mlp, _tp6);
L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr50 */


double singularhjmpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, int _tp5, int _tp6, int _tp7, int _tp8, int _tp9, int _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, int * _tp31, long _tpl31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37));

#if MLPROTOTYPES
static int _tr51( MLINK mlp)
#else
static int _tr51(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	int _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	int * _tp31;
	long _tpl31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetIntegerList( mlp, &_tp31, &_tpl31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLNewPacket(mlp) ) goto L37;

	_rp0 = singularhjmpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tpl31, _tp32, _tp33, _tp34, _tp35, _tp36, _tp37);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L37: L36: L35: L34: L33: L32: L31:	MLReleaseInteger32List( mlp, _tp31, _tpl31);
L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr51 */


double singulardoublepiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, int _tp6, int _tp7, int _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, double _tp29, double _tp30, double _tp31, int * _tp32, long _tpl32, double _tp33, double _tp34, double _tp35, double _tp36, double _tp37, double _tp38, double _tp39));

#if MLPROTOTYPES
static int _tr52( MLINK mlp)
#else
static int _tr52(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	double _tp29;
	double _tp30;
	double _tp31;
	int * _tp32;
	long _tpl32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _tp37;
	double _tp38;
	double _tp39;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetReal( mlp, &_tp29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetIntegerList( mlp, &_tp32, &_tpl32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLGetReal( mlp, &_tp37) ) goto L36;
	if ( ! MLGetReal( mlp, &_tp38) ) goto L37;
	if ( ! MLGetReal( mlp, &_tp39) ) goto L38;
	if ( ! MLNewPacket(mlp) ) goto L39;

	_rp0 = singulardoublepiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tp30, _tp31, _tp32, _tpl32, _tp33, _tp34, _tp35, _tp36, _tp37, _tp38, _tp39);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L39: L38: L37: L36: L35: L34: L33: L32:	MLReleaseInteger32List( mlp, _tp32, _tpl32);
L31: L30: L29: L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr52 */


double singularpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, int _tp6, int _tp7, int _tp8, int _tp9, int _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, int * _tp29, long _tpl29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35));

#if MLPROTOTYPES
static int _tr53( MLINK mlp)
#else
static int _tr53(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	int * _tp29;
	long _tpl29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetIntegerList( mlp, &_tp29, &_tpl29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLNewPacket(mlp) ) goto L35;

	_rp0 = singularpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tpl29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L35: L34: L33: L32: L31: L30: L29:	MLReleaseInteger32List( mlp, _tp29, _tpl29);
L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr53 */


double singulardiffpiece P(( const char * _tp1, const char * _tp2, const char * _tp3, const char * _tp4, const char * _tp5, int _tp6, int _tp7, int _tp8, int _tp9, int _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21, double _tp22, double _tp23, double _tp24, double _tp25, double _tp26, double _tp27, double _tp28, int * _tp29, long _tpl29, double _tp30, double _tp31, double _tp32, double _tp33, double _tp34, double _tp35, double _tp36));

#if MLPROTOTYPES
static int _tr54( MLINK mlp)
#else
static int _tr54(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	const char * _tp3;
	const char * _tp4;
	const char * _tp5;
	int _tp6;
	int _tp7;
	int _tp8;
	int _tp9;
	int _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _tp22;
	double _tp23;
	double _tp24;
	double _tp25;
	double _tp26;
	double _tp27;
	double _tp28;
	int * _tp29;
	long _tpl29;
	double _tp30;
	double _tp31;
	double _tp32;
	double _tp33;
	double _tp34;
	double _tp35;
	double _tp36;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetString( mlp, &_tp3) ) goto L2;
	if ( ! MLGetString( mlp, &_tp4) ) goto L3;
	if ( ! MLGetString( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLGetReal( mlp, &_tp22) ) goto L21;
	if ( ! MLGetReal( mlp, &_tp23) ) goto L22;
	if ( ! MLGetReal( mlp, &_tp24) ) goto L23;
	if ( ! MLGetReal( mlp, &_tp25) ) goto L24;
	if ( ! MLGetReal( mlp, &_tp26) ) goto L25;
	if ( ! MLGetReal( mlp, &_tp27) ) goto L26;
	if ( ! MLGetReal( mlp, &_tp28) ) goto L27;
	if ( ! MLGetIntegerList( mlp, &_tp29, &_tpl29) ) goto L28;
	if ( ! MLGetReal( mlp, &_tp30) ) goto L29;
	if ( ! MLGetReal( mlp, &_tp31) ) goto L30;
	if ( ! MLGetReal( mlp, &_tp32) ) goto L31;
	if ( ! MLGetReal( mlp, &_tp33) ) goto L32;
	if ( ! MLGetReal( mlp, &_tp34) ) goto L33;
	if ( ! MLGetReal( mlp, &_tp35) ) goto L34;
	if ( ! MLGetReal( mlp, &_tp36) ) goto L35;
	if ( ! MLNewPacket(mlp) ) goto L36;

	_rp0 = singulardiffpiece(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21, _tp22, _tp23, _tp24, _tp25, _tp26, _tp27, _tp28, _tp29, _tpl29, _tp30, _tp31, _tp32, _tp33, _tp34, _tp35, _tp36);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L36: L35: L34: L33: L32: L31: L30: L29:	MLReleaseInteger32List( mlp, _tp29, _tpl29);
L28: L27: L26: L25: L24: L23: L22: L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5:	MLReleaseString(mlp, _tp5);
L4:	MLReleaseString(mlp, _tp4);
L3:	MLReleaseString(mlp, _tp3);
L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr54 */


double diffdeltagap P(( const char * _tp1, const char * _tp2, int _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, int _tp9, int _tp10, int _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19));

#if MLPROTOTYPES
static int _tr55( MLINK mlp)
#else
static int _tr55(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	int _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	int _tp9;
	int _tp10;
	int _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetInteger( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLNewPacket(mlp) ) goto L19;

	_rp0 = diffdeltagap(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr55 */


double diffdeltagapmass P(( const char * _tp1, int _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, int _tp10, int _tp11, int _tp12, int _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19, double _tp20, double _tp21));

#if MLPROTOTYPES
static int _tr56( MLINK mlp)
#else
static int _tr56(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	int _tp10;
	int _tp11;
	int _tp12;
	int _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _tp20;
	double _tp21;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetInteger( mlp, &_tp10) ) goto L9;
	if ( ! MLGetInteger( mlp, &_tp11) ) goto L10;
	if ( ! MLGetInteger( mlp, &_tp12) ) goto L11;
	if ( ! MLGetInteger( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLGetReal( mlp, &_tp20) ) goto L19;
	if ( ! MLGetReal( mlp, &_tp21) ) goto L20;
	if ( ! MLNewPacket(mlp) ) goto L21;

	_rp0 = diffdeltagapmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19, _tp20, _tp21);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L21: L20: L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr56 */


double hyperf32exact P(( double _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr57( MLINK mlp)
#else
static int _tr57(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _rp0;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	_rp0 = hyperf32exact(_tp1, _tp2);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L2: L1: 
L0:	return res;
} /* _tr57 */


double model P(( double * _tp1, long _tpl1, double _tp2, int _tp3, double _tp4));

#if MLPROTOTYPES
static int _tr58( MLINK mlp)
#else
static int _tr58(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double * _tp1;
	long _tpl1;
	double _tp2;
	int _tp3;
	double _tp4;
	double _rp0;
	if ( ! MLGetRealList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	_rp0 = model(_tp1, _tpl1, _tp2, _tp3, _tp4);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L4: L3: L2: L1:	MLReleaseReal64List(mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr58 */


double modelunstable P(( const char * _tp1, double _tp2, double _tp3, double * _tp4, long _tpl4, double _tp5, int _tp6, int _tp7, double _tp8));

#if MLPROTOTYPES
static int _tr59( MLINK mlp)
#else
static int _tr59(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	double _tp2;
	double _tp3;
	double * _tp4;
	long _tpl4;
	double _tp5;
	int _tp6;
	int _tp7;
	double _tp8;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetRealList( mlp, &_tp4, &_tpl4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLNewPacket(mlp) ) goto L8;

	_rp0 = modelunstable(_tp1, _tp2, _tp3, _tp4, _tpl4, _tp5, _tp6, _tp7, _tp8);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L8: L7: L6: L5: L4:	MLReleaseReal64List(mlp, _tp4, _tpl4);
L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr59 */


double breitmodelunstable P(( const char * _tp1, double _tp2, double _tp3, double _tp4, double * _tp5, long _tpl5, double _tp6, int _tp7, int _tp8, double _tp9));

#if MLPROTOTYPES
static int _tr60( MLINK mlp)
#else
static int _tr60(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double * _tp5;
	long _tpl5;
	double _tp6;
	int _tp7;
	int _tp8;
	double _tp9;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetRealList( mlp, &_tp5, &_tpl5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetInteger( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLNewPacket(mlp) ) goto L9;

	_rp0 = breitmodelunstable(_tp1, _tp2, _tp3, _tp4, _tp5, _tpl5, _tp6, _tp7, _tp8, _tp9);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L9: L8: L7: L6: L5:	MLReleaseReal64List(mlp, _tp5, _tpl5);
L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr60 */


double modelunstablediff P(( const char * _tp1, double _tp2, double _tp3, double * _tp4, long _tpl4, double _tp5, int _tp6, int _tp7, double _tp8, double _tp9));

#if MLPROTOTYPES
static int _tr61( MLINK mlp)
#else
static int _tr61(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	double _tp2;
	double _tp3;
	double * _tp4;
	long _tpl4;
	double _tp5;
	int _tp6;
	int _tp7;
	double _tp8;
	double _tp9;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetRealList( mlp, &_tp4, &_tpl4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLNewPacket(mlp) ) goto L9;

	_rp0 = modelunstablediff(_tp1, _tp2, _tp3, _tp4, _tpl4, _tp5, _tp6, _tp7, _tp8, _tp9);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L9: L8: L7: L6: L5: L4:	MLReleaseReal64List(mlp, _tp4, _tpl4);
L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr61 */


double modeldiff P(( double * _tp1, long _tpl1, double _tp2, int _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr62( MLINK mlp)
#else
static int _tr62(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double * _tp1;
	long _tpl1;
	double _tp2;
	int _tp3;
	double _tp4;
	double _tp5;
	double _rp0;
	if ( ! MLGetRealList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	_rp0 = modeldiff(_tp1, _tpl1, _tp2, _tp3, _tp4, _tp5);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L5: L4: L3: L2: L1:	MLReleaseReal64List(mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr62 */


double breitmodel P(( double * _tp1, long _tpl1, double _tp2, double _tp3, int _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr63( MLINK mlp)
#else
static int _tr63(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double * _tp1;
	long _tpl1;
	double _tp2;
	double _tp3;
	int _tp4;
	double _tp5;
	double _rp0;
	if ( ! MLGetRealList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	_rp0 = breitmodel(_tp1, _tpl1, _tp2, _tp3, _tp4, _tp5);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L5: L4: L3: L2: L1:	MLReleaseReal64List(mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr63 */


double breitmodeldiff P(( double * _tp1, long _tpl1, double _tp2, double _tp3, int _tp4, double _tp5, double _tp6));

#if MLPROTOTYPES
static int _tr64( MLINK mlp)
#else
static int _tr64(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double * _tp1;
	long _tpl1;
	double _tp2;
	double _tp3;
	int _tp4;
	double _tp5;
	double _tp6;
	double _rp0;
	if ( ! MLGetRealList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLNewPacket(mlp) ) goto L6;

	_rp0 = breitmodeldiff(_tp1, _tpl1, _tp2, _tp3, _tp4, _tp5, _tp6);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L6: L5: L4: L3: L2: L1:	MLReleaseReal64List(mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr64 */


double taylor P(( double * _tp1, long _tpl1, double _tp2, int _tp3));

#if MLPROTOTYPES
static int _tr65( MLINK mlp)
#else
static int _tr65(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double * _tp1;
	long _tpl1;
	double _tp2;
	int _tp3;
	double _rp0;
	if ( ! MLGetRealList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	_rp0 = taylor(_tp1, _tpl1, _tp2, _tp3);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L3: L2: L1:	MLReleaseReal64List(mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr65 */


double momentmodel P(( double * _tp1, long _tpl1, double _tp2, int _tp3));

#if MLPROTOTYPES
static int _tr66( MLINK mlp)
#else
static int _tr66(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double * _tp1;
	long _tpl1;
	double _tp2;
	int _tp3;
	double _rp0;
	if ( ! MLGetRealList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	_rp0 = momentmodel(_tp1, _tpl1, _tp2, _tp3);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L3: L2: L1:	MLReleaseReal64List(mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr66 */


double modelpiece P(( int * _tp1, long _tpl1, double _tp2, int _tp3, double _tp4));

#if MLPROTOTYPES
static int _tr67( MLINK mlp)
#else
static int _tr67(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int * _tp1;
	long _tpl1;
	double _tp2;
	int _tp3;
	double _tp4;
	double _rp0;
	if ( ! MLGetIntegerList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	_rp0 = modelpiece(_tp1, _tpl1, _tp2, _tp3, _tp4);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L4: L3: L2: L1:	MLReleaseInteger32List( mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr67 */


double taylorpiece P(( int * _tp1, long _tpl1, double _tp2, int _tp3));

#if MLPROTOTYPES
static int _tr68( MLINK mlp)
#else
static int _tr68(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int * _tp1;
	long _tpl1;
	double _tp2;
	int _tp3;
	double _rp0;
	if ( ! MLGetIntegerList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	_rp0 = taylorpiece(_tp1, _tpl1, _tp2, _tp3);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L3: L2: L1:	MLReleaseInteger32List( mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr68 */


double hyper2f1 P(( double _tp1, double _tp2, double _tp3, double _tp4));

#if MLPROTOTYPES
static int _tr69( MLINK mlp)
#else
static int _tr69(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _rp0;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	_rp0 = hyper2f1(_tp1, _tp2, _tp3, _tp4);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L4: L3: L2: L1: 
L0:	return res;
} /* _tr69 */


double intecorre P(( double _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr70( MLINK mlp)
#else
static int _tr70(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _rp0;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	_rp0 = intecorre(_tp1, _tp2, _tp3);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L3: L2: L1: 
L0:	return res;
} /* _tr70 */


void anomdim P(( const char * _tp1, int _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr71( MLINK mlp)
#else
static int _tr71(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	double _tp3;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	anomdim(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr71 */


void scoef P(( const char * _tp1, int _tp2));

#if MLPROTOTYPES
static int _tr72( MLINK mlp)
#else
static int _tr72(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	scoef(_tp1, _tp2);

	res = 1;
L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr72 */


double polylog P(( int _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr73( MLINK mlp)
#else
static int _tr73(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	_rp0 = polylog(_tp1, _tp2);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L2: L1: 
L0:	return res;
} /* _tr73 */


double dilog P(( double _tp1));

#if MLPROTOTYPES
static int _tr74( MLINK mlp)
#else
static int _tr74(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _rp0;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	_rp0 = dilog(_tp1);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L1: 
L0:	return res;
} /* _tr74 */


double pfq P(( double * _tp1, long _tpl1, double * _tp2, long _tpl2, double _tp3));

#if MLPROTOTYPES
static int _tr75( MLINK mlp)
#else
static int _tr75(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double * _tp1;
	long _tpl1;
	double * _tp2;
	long _tpl2;
	double _tp3;
	double _rp0;
	if ( ! MLGetRealList( mlp, &_tp1, &_tpl1) ) goto L0;
	if ( ! MLGetRealList( mlp, &_tp2, &_tpl2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	_rp0 = pfq(_tp1, _tpl1, _tp2, _tpl2, _tp3);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L3: L2:	MLReleaseReal64List(mlp, _tp2, _tpl2);
L1:	MLReleaseReal64List(mlp, _tp1, _tpl1);

L0:	return res;
} /* _tr75 */


double elliptic3 P(( double _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr76( MLINK mlp)
#else
static int _tr76(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	double _tp3;
	double _rp0;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	_rp0 = elliptic3(_tp1, _tp2, _tp3);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L3: L2: L1: 
L0:	return res;
} /* _tr76 */


double nglfunction P(( int _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr77( MLINK mlp)
#else
static int _tr77(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	_rp0 = nglfunction(_tp1, _tp2);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L2: L1: 
L0:	return res;
} /* _tr77 */


void nglsoft P(( int _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr78( MLINK mlp)
#else
static int _tr78(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	nglsoft(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr78 */


void complexpolylog P(( int _tp1, double _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr79( MLINK mlp)
#else
static int _tr79(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	complexpolylog(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1: 
L0:	return res;
} /* _tr79 */


void cli2 P(( double _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr80( MLINK mlp)
#else
static int _tr80(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	cli2(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr80 */


void thrustns1loop P(( double _tp1));

#if MLPROTOTYPES
static int _tr81( MLINK mlp)
#else
static int _tr81(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLNewPacket(mlp) ) goto L1;

	thrustns1loop(_tp1);

	res = 1;
L1: 
L0:	return res;
} /* _tr81 */


double fomass P(( const char * _tp1, const char * _tp2, double _tp3, double _tp4, double _tp5, double _tp6, double _tp7, double _tp8));

#if MLPROTOTYPES
static int _tr82( MLINK mlp)
#else
static int _tr82(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLNewPacket(mlp) ) goto L8;

	_rp0 = fomass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L8: L7: L6: L5: L4: L3: L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr82 */


void thrustns2loop P(( double _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr83( MLINK mlp)
#else
static int _tr83(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	thrustns2loop(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr83 */


void cli3 P(( double _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr84( MLINK mlp)
#else
static int _tr84(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	double _tp1;
	double _tp2;
	if ( ! MLGetReal( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	cli3(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr84 */


void delta P(( const char * _tp1, int _tp2, double _tp3, double _tp4));

#if MLPROTOTYPES
static int _tr85( MLINK mlp)
#else
static int _tr85(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	double _tp3;
	double _tp4;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	delta(_tp1, _tp2, _tp3, _tp4);

	res = 1;
L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr85 */


void deltagap P(( const char * _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15));

#if MLPROTOTYPES
static int _tr86( MLINK mlp)
#else
static int _tr86(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLNewPacket(mlp) ) goto L15;

	deltagap(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15);

	res = 1;
L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr86 */


void coefmat P(( const char * _tp1, int _tp2, double _tp3));

#if MLPROTOTYPES
static int _tr87( MLINK mlp)
#else
static int _tr87(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	double _tp3;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLNewPacket(mlp) ) goto L3;

	coefmat(_tp1, _tp2, _tp3);

	res = 1;
L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr87 */


double wtilde P(( int _tp1, int _tp2, double * _tp3, long _tpl3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr88( MLINK mlp)
#else
static int _tr88(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double * _tp3;
	long _tpl3;
	double _tp4;
	double _tp5;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetRealList( mlp, &_tp3, &_tpl3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	_rp0 = wtilde(_tp1, _tp2, _tp3, _tpl3, _tp4, _tp5);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L5: L4: L3:	MLReleaseReal64List(mlp, _tp3, _tpl3);
L2: L1: 
L0:	return res;
} /* _tr88 */


double ktilde P(( int _tp1, int _tp2, double * _tp3, long _tpl3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr89( MLINK mlp)
#else
static int _tr89(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double * _tp3;
	long _tpl3;
	double _tp4;
	double _tp5;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetRealList( mlp, &_tp3, &_tpl3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	_rp0 = ktilde(_tp1, _tp2, _tp3, _tpl3, _tp4, _tp5);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L5: L4: L3:	MLReleaseReal64List(mlp, _tp3, _tpl3);
L2: L1: 
L0:	return res;
} /* _tr89 */


double alphaqcd P(( const char * _tp1, const char * _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14));

#if MLPROTOTYPES
static int _tr90( MLINK mlp)
#else
static int _tr90(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLNewPacket(mlp) ) goto L14;

	_rp0 = alphaqcd(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr90 */


void alphacomplex P(( const char * _tp1, const char * _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15));

#if MLPROTOTYPES
static int _tr91( MLINK mlp)
#else
static int _tr91(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLNewPacket(mlp) ) goto L15;

	alphacomplex(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15);

	res = 1;
L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr91 */


double msbarmass P(( int _tp1, int _tp2, int _tp3, int _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13));

#if MLPROTOTYPES
static int _tr92( MLINK mlp)
#else
static int _tr92(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLNewPacket(mlp) ) goto L13;

	_rp0 = msbarmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr92 */


double polemass P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14));

#if MLPROTOTYPES
static int _tr93( MLINK mlp)
#else
static int _tr93(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLNewPacket(mlp) ) goto L14;

	_rp0 = polemass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr93 */


double msbarmasslow P(( int _tp1, int _tp2, int _tp3, int _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13));

#if MLPROTOTYPES
static int _tr94( MLINK mlp)
#else
static int _tr94(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLNewPacket(mlp) ) goto L13;

	_rp0 = msbarmasslow(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr94 */


double msrmass P(( int _tp1, int _tp2, int _tp3, int _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14));

#if MLPROTOTYPES
static int _tr95( MLINK mlp)
#else
static int _tr95(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLNewPacket(mlp) ) goto L14;

	_rp0 = msrmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr95 */


double mmfrommsr P(( int _tp1, int _tp2, int _tp3, int _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14));

#if MLPROTOTYPES
static int _tr96( MLINK mlp)
#else
static int _tr96(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLNewPacket(mlp) ) goto L14;

	_rp0 = mmfrommsr(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr96 */


double mmfrommsrnatural P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15));

#if MLPROTOTYPES
static int _tr97( MLINK mlp)
#else
static int _tr97(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLNewPacket(mlp) ) goto L15;

	_rp0 = mmfrommsrnatural(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr97 */


double msrnaturalmass P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15));

#if MLPROTOTYPES
static int _tr98( MLINK mlp)
#else
static int _tr98(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLNewPacket(mlp) ) goto L15;

	_rp0 = msrnaturalmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr98 */


double jetmass P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16));

#if MLPROTOTYPES
static int _tr99( MLINK mlp)
#else
static int _tr99(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLNewPacket(mlp) ) goto L16;

	_rp0 = jetmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr99 */


double mmfromjetmass P(( int _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16));

#if MLPROTOTYPES
static int _tr100( MLINK mlp)
#else
static int _tr100(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLNewPacket(mlp) ) goto L16;

	_rp0 = mmfromjetmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr100 */


void deltamsbar P(( int _tp1, int _tp2, int _tp3, int _tp4, double _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13));

#if MLPROTOTYPES
static int _tr101( MLINK mlp)
#else
static int _tr101(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	double _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLNewPacket(mlp) ) goto L13;

	deltamsbar(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13);

	res = 1;
L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr101 */


double rhad P(( const char * _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15));

#if MLPROTOTYPES
static int _tr102( MLINK mlp)
#else
static int _tr102(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLNewPacket(mlp) ) goto L15;

	_rp0 = rhad(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr102 */


double rhadmass P(( const char * _tp1, const char * _tp2, int _tp3, int _tp4, int _tp5, int _tp6, int _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14, double _tp15, double _tp16, double _tp17, double _tp18, double _tp19));

#if MLPROTOTYPES
static int _tr103( MLINK mlp)
#else
static int _tr103(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	const char * _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	int _tp6;
	int _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _tp15;
	double _tp16;
	double _tp17;
	double _tp18;
	double _tp19;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetString( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetInteger( mlp, &_tp6) ) goto L5;
	if ( ! MLGetInteger( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLGetReal( mlp, &_tp15) ) goto L14;
	if ( ! MLGetReal( mlp, &_tp16) ) goto L15;
	if ( ! MLGetReal( mlp, &_tp17) ) goto L16;
	if ( ! MLGetReal( mlp, &_tp18) ) goto L17;
	if ( ! MLGetReal( mlp, &_tp19) ) goto L18;
	if ( ! MLNewPacket(mlp) ) goto L19;

	_rp0 = rhadmass(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14, _tp15, _tp16, _tp17, _tp18, _tp19);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L19: L18: L17: L16: L15: L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2:	MLReleaseString(mlp, _tp2);
L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr103 */


double lambdaqcd P(( const char * _tp1, int _tp2, int _tp3, int _tp4, int _tp5, double _tp6, double _tp7, double _tp8, double _tp9, double _tp10, double _tp11, double _tp12, double _tp13, double _tp14));

#if MLPROTOTYPES
static int _tr104( MLINK mlp)
#else
static int _tr104(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	const char * _tp1;
	int _tp2;
	int _tp3;
	int _tp4;
	int _tp5;
	double _tp6;
	double _tp7;
	double _tp8;
	double _tp9;
	double _tp10;
	double _tp11;
	double _tp12;
	double _tp13;
	double _tp14;
	double _rp0;
	if ( ! MLGetString( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetInteger( mlp, &_tp4) ) goto L3;
	if ( ! MLGetInteger( mlp, &_tp5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetReal( mlp, &_tp7) ) goto L6;
	if ( ! MLGetReal( mlp, &_tp8) ) goto L7;
	if ( ! MLGetReal( mlp, &_tp9) ) goto L8;
	if ( ! MLGetReal( mlp, &_tp10) ) goto L9;
	if ( ! MLGetReal( mlp, &_tp11) ) goto L10;
	if ( ! MLGetReal( mlp, &_tp12) ) goto L11;
	if ( ! MLGetReal( mlp, &_tp13) ) goto L12;
	if ( ! MLGetReal( mlp, &_tp14) ) goto L13;
	if ( ! MLNewPacket(mlp) ) goto L14;

	_rp0 = lambdaqcd(_tp1, _tp2, _tp3, _tp4, _tp5, _tp6, _tp7, _tp8, _tp9, _tp10, _tp11, _tp12, _tp13, _tp14);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L14: L13: L12: L11: L10: L9: L8: L7: L6: L5: L4: L3: L2: L1:	MLReleaseString(mlp, _tp1);

L0:	return res;
} /* _tr104 */


void kernels P(( int _tp1, double _tp2, double _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr105( MLINK mlp)
#else
static int _tr105(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	kernels(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = 1;
L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr105 */


void gammaderlist P(( int _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr106( MLINK mlp)
#else
static int _tr106(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	gammaderlist(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr106 */


void polygamma P(( int _tp1, double _tp2));

#if MLPROTOTYPES
static int _tr107( MLINK mlp)
#else
static int _tr107(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	double _tp2;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetReal( mlp, &_tp2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	polygamma(_tp1, _tp2);

	res = 1;
L2: L1: 
L0:	return res;
} /* _tr107 */


void nglkernels P(( int _tp1, int _tp2, int _tp3, double _tp4, double * _tp5, long _tpl5, double _tp6, double * _tp7, long _tpl7));

#if MLPROTOTYPES
static int _tr108( MLINK mlp)
#else
static int _tr108(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	int _tp3;
	double _tp4;
	double * _tp5;
	long _tpl5;
	double _tp6;
	double * _tp7;
	long _tpl7;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetInteger( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetRealList( mlp, &_tp5, &_tpl5) ) goto L4;
	if ( ! MLGetReal( mlp, &_tp6) ) goto L5;
	if ( ! MLGetRealList( mlp, &_tp7, &_tpl7) ) goto L6;
	if ( ! MLNewPacket(mlp) ) goto L7;

	nglkernels(_tp1, _tp2, _tp3, _tp4, _tp5, _tpl5, _tp6, _tp7, _tpl7);

	res = 1;
L7:	MLReleaseReal64List(mlp, _tp7, _tpl7);
L6: L5:	MLReleaseReal64List(mlp, _tp5, _tpl5);
L4: L3: L2: L1: 
L0:	return res;
} /* _tr108 */


double nglintegral P(( int _tp1, int _tp2, double _tp3, double _tp4));

#if MLPROTOTYPES
static int _tr109( MLINK mlp)
#else
static int _tr109(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double _tp3;
	double _tp4;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLNewPacket(mlp) ) goto L4;

	_rp0 = nglintegral(_tp1, _tp2, _tp3, _tp4);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L4: L3: L2: L1: 
L0:	return res;
} /* _tr109 */


double ngldoubleintegral P(( int _tp1, int _tp2, double _tp3, double _tp4, double _tp5));

#if MLPROTOTYPES
static int _tr110( MLINK mlp)
#else
static int _tr110(mlp) MLINK mlp;
#endif
{
	int	res = 0;
	int _tp1;
	int _tp2;
	double _tp3;
	double _tp4;
	double _tp5;
	double _rp0;
	if ( ! MLGetInteger( mlp, &_tp1) ) goto L0;
	if ( ! MLGetInteger( mlp, &_tp2) ) goto L1;
	if ( ! MLGetReal( mlp, &_tp3) ) goto L2;
	if ( ! MLGetReal( mlp, &_tp4) ) goto L3;
	if ( ! MLGetReal( mlp, &_tp5) ) goto L4;
	if ( ! MLNewPacket(mlp) ) goto L5;

	_rp0 = ngldoubleintegral(_tp1, _tp2, _tp3, _tp4, _tp5);

	res = MLAbort ?
		MLPutFunction( mlp, "Abort", 0) : MLPutReal( mlp, _rp0);
L5: L4: L3: L2: L1: 
L0:	return res;
} /* _tr110 */


static struct func {
	int   f_nargs;
	int   manual;
	int   (*f_func)P((MLINK));
	const char  *f_name;
	} _tramps[111] = {
		{ 5, 0, _tr0, "ewfactors" },
		{ 3, 0, _tr1, "legendrelist" },
		{ 2, 0, _tr2, "qlegendrelist" },
		{ 2, 0, _tr3, "gammar" },
		{46, 0, _tr4, "masslessprof" },
		{33, 0, _tr5, "findorigin" },
		{45, 0, _tr6, "masslessprofpiece" },
		{45, 0, _tr7, "masslessprofpiecelist" },
		{46, 0, _tr8, "masslesspiecebin" },
		{46, 0, _tr9, "masslessprofdiffpiece" },
		{46, 0, _tr10, "masslessproflist" },
		{47, 0, _tr11, "masslessbinlist" },
		{47, 0, _tr12, "masslessprofdiff" },
		{47, 0, _tr13, "masslessmoment" },
		{16, 0, _tr14, "profiles" },
		{22, 0, _tr15, "profilesmass" },
		{ 6, 0, _tr16, "mctop" },
		{ 7, 0, _tr17, "breitunstable" },
		{ 3, 0, _tr18, "deltamctop" },
		{19, 0, _tr19, "massintersection" },
		{40, 0, _tr20, "nsmasspiece" },
		{41, 0, _tr21, "nsmassdiffpiece" },
		{41, 0, _tr22, "nsmass" },
		{42, 0, _tr23, "nsmassdiff" },
		{41, 0, _tr24, "hjmnsmass" },
		{42, 0, _tr25, "hjmnsmassdiff" },
		{36, 0, _tr26, "singular" },
		{35, 0, _tr27, "singularlist" },
		{37, 0, _tr28, "singulardiff" },
		{62, 0, _tr29, "massiveprof" },
		{43, 0, _tr30, "massorigin" },
		{62, 0, _tr31, "massiveprofpiece" },
		{62, 0, _tr32, "massiveprofpiecelist" },
		{63, 0, _tr33, "massivepiecebin" },
		{63, 0, _tr34, "massiveprofdiffpiece" },
		{63, 0, _tr35, "massiveprofdiff" },
		{62, 0, _tr36, "massiveproflist" },
		{63, 0, _tr37, "massivebinlist" },
		{63, 0, _tr38, "massivemoment" },
		{50, 0, _tr39, "singularmass" },
		{51, 0, _tr40, "singularmassdiff" },
		{41, 0, _tr41, "massnondist" },
		{42, 0, _tr42, "massnondistdiff" },
		{40, 0, _tr43, "massnondistpiece" },
		{41, 0, _tr44, "massnondistdiffpiece" },
		{50, 0, _tr45, "singularmasspiece" },
		{51, 0, _tr46, "singularmassdiffpiece" },
		{39, 0, _tr47, "singularhjm" },
		{36, 0, _tr48, "singularhjm1d" },
		{36, 0, _tr49, "singularhjm1dpiece" },
		{41, 0, _tr50, "singulardouble" },
		{37, 0, _tr51, "singularhjmpiece" },
		{39, 0, _tr52, "singulardoublepiece" },
		{35, 0, _tr53, "singularpiece" },
		{36, 0, _tr54, "singulardiffpiece" },
		{19, 0, _tr55, "diffdeltagap" },
		{21, 0, _tr56, "diffdeltagapmass" },
		{ 2, 0, _tr57, "hyperf32exact" },
		{ 4, 0, _tr58, "model" },
		{ 8, 0, _tr59, "modelunstable" },
		{ 9, 0, _tr60, "breitmodelunstable" },
		{ 9, 0, _tr61, "modelunstablediff" },
		{ 5, 0, _tr62, "modeldiff" },
		{ 5, 0, _tr63, "breitmodel" },
		{ 6, 0, _tr64, "breitmodeldiff" },
		{ 3, 0, _tr65, "taylor" },
		{ 3, 0, _tr66, "momentmodel" },
		{ 4, 0, _tr67, "modelpiece" },
		{ 3, 0, _tr68, "taylorpiece" },
		{ 4, 0, _tr69, "hyper2f1" },
		{ 3, 0, _tr70, "intecorre" },
		{ 3, 0, _tr71, "anomdim" },
		{ 2, 0, _tr72, "scoef" },
		{ 2, 0, _tr73, "polylog" },
		{ 1, 0, _tr74, "dilog" },
		{ 3, 0, _tr75, "pfq" },
		{ 3, 0, _tr76, "elliptic3" },
		{ 2, 0, _tr77, "nglfunction" },
		{ 3, 0, _tr78, "nglsoft" },
		{ 3, 0, _tr79, "complexpolylog" },
		{ 2, 0, _tr80, "cli2" },
		{ 1, 0, _tr81, "thrustns1loop" },
		{ 8, 0, _tr82, "fomass" },
		{ 2, 0, _tr83, "thrustns2loop" },
		{ 2, 0, _tr84, "cli3" },
		{ 4, 0, _tr85, "delta" },
		{15, 0, _tr86, "deltagap" },
		{ 3, 0, _tr87, "coefmat" },
		{ 5, 0, _tr88, "wtilde" },
		{ 5, 0, _tr89, "ktilde" },
		{14, 0, _tr90, "alphaqcd" },
		{15, 0, _tr91, "alphacomplex" },
		{13, 0, _tr92, "msbarmass" },
		{14, 0, _tr93, "polemass" },
		{13, 0, _tr94, "msbarmasslow" },
		{14, 0, _tr95, "msrmass" },
		{14, 0, _tr96, "mmfrommsr" },
		{15, 0, _tr97, "mmfrommsrnatural" },
		{15, 0, _tr98, "msrnaturalmass" },
		{16, 0, _tr99, "jetmass" },
		{16, 0, _tr100, "mmfromjetmass" },
		{13, 0, _tr101, "deltamsbar" },
		{15, 0, _tr102, "rhad" },
		{19, 0, _tr103, "rhadmass" },
		{14, 0, _tr104, "lambdaqcd" },
		{ 5, 0, _tr105, "kernels" },
		{ 2, 0, _tr106, "gammaderlist" },
		{ 2, 0, _tr107, "polygamma" },
		{ 7, 0, _tr108, "nglkernels" },
		{ 4, 0, _tr109, "nglintegral" },
		{ 5, 0, _tr110, "ngldoubleintegral" }
		};

static const char* evalstrs[] = {
	"BeginPackage[\"Caliper`\"]",
	(const char*)0,
	"Print[\"     Package for Massive and Massless Event Shapes \"]",
	(const char*)0,
	"Print[\"     Author:            Vicent Mateu               \"]",
	(const char*)0,
	"Print[\"     Last modification: 23 - 02 - 2017             \"]",
	(const char*)0,
	"Print[\"     Version:           test 1                     \"]",
	(const char*)0,
	"mZdef            = 91.187",
	(const char*)0,
	"Gammadef         = 1553.0647546066",
	(const char*)0,
	"gammaZdef        = 2.4952",
	(const char*)0,
	"sin2ThetaWdef    = 0.23119",
	(const char*)0,
	"gammaZPythia     = 2.5042",
	(const char*)0,
	"sin2ThetaWPythia = 0.2312",
	(const char*)0,
	"LegendreList::usage = \"LegendreList[n, k, x] computes the of the",
	" first n + 1 k-th derivative of the Legendre Polynomials\"",
	(const char*)0,
	"QLegendreList::usage = \"QLegendreList[n, x] computes the of the ",
	"first n + 1 Legendre Polynomial\"",
	(const char*)0,
	"BreitUnstable::usage = \"BreitUnstable[shape, mt, Q, gamma, n, k,",
	" x] computes the LO distribution for unstable tops convoluted wi",
	"th a BreitWigner\"",
	(const char*)0,
	"MCtop::usage = \"MCtop[shape, mt, Q, n, k, x] computes the LO dis",
	"tribution for unstable tops\"",
	(const char*)0,
	"DeltaMCtop::usage = \"DeltaMCtop[shape, mt, Q] computes the LO di",
	"stribution for unstable tops\"",
	(const char*)0,
	"pFq::usage = \"pFq[a,b,z] computes the generalized hypergeometric",
	" function\"",
	(const char*)0,
	"FindOrigin::usage = \"FindOrigin[shape, gap, orderAlpha, runAlpha",
	", order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, ",
	"Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eR, R0, mu",
	"R0, delta0, h], finds the origin for massless Event Shapes\"",
	(const char*)0,
	"MassOrigin::usage = \"MassOrigin[shape, EShape, gap, scheme, orde",
	"rAlpha, runAlpha, orderMass, runMass, order, run, nf, mZ, amZ, m",
	"T, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0, de",
	"ltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH,",
	" eS, eJ, mass, muM, R0, muR0, del0, h], finds the origin of Mass",
	"ive Event Shapes\"",
	(const char*)0,
	"MasslessPieceBin::usage = \"MasslessPieceBin[terms, hard, shape, ",
	"setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3",
	", s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, ",
	"Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, ",
	"lambda, R0, muR0, delta0, h, tauList], computes bins including p",
	"rofiles for massless cross section\"",
	(const char*)0,
	"MasslessBinList::usage = \"MasslessBinList[terms, hard, shape, se",
	"tup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, ",
	"s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Ra",
	"t0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambd",
	"a, R0, muR0, delta0, h, tauList], computes bins including profil",
	"es for massless cross section\"",
	(const char*)0,
	"MasslessProfList::usage = \"MasslessProfList[terms, hard, shape, ",
	"setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3",
	", s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, ",
	"Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lam",
	"bda, R0, muR0, delta0, h, tauList], computes the cross section i",
	"ncluding profiles for massless cross section\"",
	(const char*)0,
	"MassiveBinList::usage = \"MassiveBinList[terms, hard, shape, ESha",
	"pe, setup, gap, space, cum, scheme, abs, current, xi, xiB, order",
	"Alpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3,",
	" mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, be",
	"ta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slop",
	"e, cnt, eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0, m",
	"uR0, del0, h, gammaZ, sin2ThetaW, tauList], computes bins includ",
	"ing profiles for massive cross section\"",
	(const char*)0,
	"MassiveProfList::usage = \"MassiveProfList[terms, hard, shape, ES",
	"hape, setup, gap, space, cum, scheme, abs, current, xi, xiB, ord",
	"erAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G",
	"3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, ",
	"beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, sl",
	"ope, cnt, eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0,",
	" muR0, del0, h, gammaZ, sin2ThetaW, tauList], computes the cross",
	" section including profiles for massive cross section\"",
	(const char*)0,
	"MassiveProf::usage = \"MassiveProf[terms, hard, shape, EShape, se",
	"tup, gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha,",
	" runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, a",
	"mZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu",
	"0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt",
	", eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0, muR0, d",
	"el0, h, gammaZ, sin2ThetaW, tau], computes the cross section inc",
	"luding profiles for massive cross section\"",
	(const char*)0,
	"MassiveProfPiece::usage = \"MassiveProfPiece[terms, hard, shape, ",
	"EShape, setup, gap, space, cum, scheme, abs, current, xi, xiB, o",
	"rderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3,",
	" G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q",
	", beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, ",
	"slope, cnt, eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, ",
	"muR0, del0, h, gammaZ, sin2ThetaW, tau], computes the cross sect",
	"ion including profiles for massive cross section\"",
	(const char*)0,
	"MassiveProfPieceList::usage = \"MassiveProfPieceList[terms, hard,",
	" shape, EShape, setup, gap, space, cum, scheme, abs, current, xi",
	", xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf,",
	" j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLa",
	"mbda2, Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, ",
	"t2, ts, slope, cnt, eH, eS, eJ, mass, muM, ns, width, clen, lamb",
	"da, R0, muR0, del0, h, gammaZ, sin2ThetaW, tau], computes the cr",
	"oss section including profiles for massive cross section\"",
	(const char*)0,
	"MassivePieceBin::usage = \"MassivePieceBin[terms, hard, shape, ES",
	"hape, setup, gap, space, cum, scheme, abs, current, xi, xiB, ord",
	"erAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G",
	"3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, ",
	"beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, sl",
	"ope, cnt, eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, mu",
	"R0, del0, h, gammaZ, sin2ThetaW, tauList], computes the cross se",
	"ction including profiles for massive cross section\"",
	(const char*)0,
	"MassiveMoment::usage = \"MassiveMoment[terms, hard, shape, EShape",
	", setup, gap, space, scheme, abs, current, xi, xiB, orderAlpha, ",
	"runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, am",
	"Z, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, beta, mu0",
	", deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,",
	" eH, eS, eJ, mass, muM, ns, width, c, clen, lambda, R0, muR0, de",
	"l0, h, gammaZ, sin2ThetaW, tau, tau2, pow], computes the moment ",
	"including profiles for massive cross section\"",
	(const char*)0,
	"MasslessMoment::usage = \"MasslessMoment[terms, hard, shape, setu",
	"p, gap, space, orderAlpha, runAlpha, order, run, nf, j3, s3, G3,",
	" mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0,",
	" n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0, ",
	"muR0, delta0, h, tau, tau2, pow], computes moments including pro",
	"files for massless cross section\"",
	(const char*)0,
	"MasslessProf::usage = \"MasslessProf[terms, hard, shape, setup, g",
	"ap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3",
	", mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0, n0",
	", n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0,",
	" muR0, delta0, h, tau], computes the cross section including pro",
	"files for massless cross section\"",
	(const char*)0,
	"MasslessProfPiece::usage = \"MasslessProfPiece[terms, hard, shape",
	", gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3,",
	" G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0, Rat0,",
	" n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen, lambd",
	"a, R0, muR0, delta0, h, tau], computes the cross section includi",
	"ng profiles for massless cross section\"",
	(const char*)0,
	"MasslessProfPieceList::usage = \"MasslessProfPieceList[terms, har",
	"d, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf,",
	" j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu",
	"0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, cle",
	"n, lambda, R0, muR0, delta0, h, tau], computes the cross section",
	" including profiles for massless cross section\"",
	(const char*)0,
	"MassIntersection::usage = \"MassIntersection[Q, beta, mu0, delLam",
	"b, n0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass,",
	" muM, def, EShape], computes the Renormalization scales for mass",
	"ive event shapes\"",
	(const char*)0,
	"ProfilesMass::usage = \"ProfilesMass[Q, beta, mu0, delLamb, R0, n",
	"0, delta0, n1, delta1, t2, ts, slope, cnt, eH, eS, eJ, mass, muM",
	", ns, def, EShape, tau], computes the Renormalization scales for",
	" massive event shapes\"",
	(const char*)0,
	"Profiles::usage = \"Profiles[Q, mu0, R0, n0, n1, t2, tR, ts, slop",
	"e, cnt, eH, eS, eJ, eR, ns, tau], computes the Renormalization s",
	"cales for massless event shapes\"",
	(const char*)0,
	"GammaR::usage = \"GammaR[str, nf] computes the soft R anomalous D",
	"imension\"",
	(const char*)0,
	"DiLog::usage = \"DiLog[n,z] computes the polylogarithm\"",
	(const char*)0,
	"Elliptic3::usage = \"Elliptic3[psi, k, c] computes the elliptic f",
	"unction of the third kind, using Carlson forms\"",
	(const char*)0,
	"Polylog::usage = \"Polylog[n,z] computes the polylogarithm\"",
	(const char*)0,
	"NGLFunction::usage = \"NGLFunction[n,z] computes the non-global p",
	"art of the soft function\"",
	(const char*)0,
	"ComplexPolylog::usage = \"ComplexPolylog[n,z] computes the polylo",
	"garithm\"",
	(const char*)0,
	"NGLSoft::usage = \"NGLSoft[nf,z] computes the NGL function in com",
	"plex space\"",
	(const char*)0,
	"CLi2::usage = \"Cli2[z] computes the complex dilogarithm\"",
	(const char*)0,
	"CLi3::usage = \"Cli3[z] computes the complex trilogarithm\"",
	(const char*)0,
	"Scoef::usage = \"Scoef[str, nf] computes the soft R anomalous Dim",
	"ension\"",
	(const char*)0,
	"AnomDim::usage = \"AnomDim[str, nf, G4] computes the QCD anomalou",
	"s function\"",
	(const char*)0,
	"Delta::usage = \"Delta[str, nf, mu, R] computes the soft renormal",
	"on subtractions\"",
	(const char*)0,
	"DeltaGap::usage = \"DeltaGap[str, orderAlpha, runAlpha, runMass, ",
	"nf, mZ, aMz, mT, muT, mB, muB, mC, muC, mu, R] computes the soft",
	" renormalon subtractions\"",
	(const char*)0,
	"CoefMat::usage = \"CoefMat[str, nf, s3] computes the hard, soft a",
	"nd jet matrix elements\"",
	(const char*)0,
	"wTilde::usage = \"wTilde[order, nf, gamma, a0, a1] computes wTild",
	"e for a given anomalous dimension gamma\"",
	(const char*)0,
	"kTilde::usage = \"kTilde[order, nf, gamma, a0, a1] computes kTild",
	"e for a given anomalous dimension gamma\"",
	(const char*)0,
	"AlphaQCD::usage = \"AlphaQCD[scheme, method, order, run, nf, Mz, ",
	"aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of alph",
	"a with flavor matching.\"",
	(const char*)0,
	"AlphaComplex::usage = \"AlphaComplex[scheme, method, order, run, ",
	"nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running",
	" of alpha with flavor matching.\"",
	(const char*)0,
	"MSbarMass::usage = \"MSbarMass[order, runAlpha, run, nf, Mz, aMz,",
	" mT, muT, mB, muB, mC, muC, mu] computes the running of the quar",
	"k masses with flavor matching.\"",
	(const char*)0,
	"PoleMass::usage = \"PoleMass[orderAlpha, runAlpha, order, run, nf",
	", Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running o",
	"f the quark masses with flavor matching.\"",
	(const char*)0,
	"MSbarMassLow::usage = \"MSbarMassLow[order, runAlpha, run, nf, Mz",
	", aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of th",
	"e quark masses with flavor matching below the mass.\"",
	(const char*)0,
	"MSRMass::usage = \"MSRMass[order, runAlpha, run, nf, Mz, aMz, mT,",
	" muT, mB, muB, mC, muC, muLambda, R] computes the MSR practical ",
	"definition running of the quark masses with flavor matching.\"",
	(const char*)0,
	"mmfromMSR::usage = \"mmfromMSR[order, runAlpha, run, nf, Mz, aMz,",
	" mT, muT, mB, muB, mC, muC, muLambda, R] computes the MSR practi",
	"cal definition running of the quark masses with flavor matching.",
	"\"",
	(const char*)0,
	"mmfromMSRNatural::usage = \"mmfromMSRNatural[orderAlpha, runAlpha",
	", order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, ",
	"R] computes the MSR natural definition running of the quark mass",
	"es with flavor matching.\"",
	(const char*)0,
	"MSRNaturalMass::usage = \"MSRNaturalMass[orderAlpha, runAlpha, or",
	"der, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, R] c",
	"omputes the MSR natural definition running of the quark masses w",
	"ith flavor matching.\"",
	(const char*)0,
	"Rhad::usage = \"Rhad[scheme, orderAlpha, runAlpha, order, nf, Mz,",
	" aMz, mT, muT, mB, muB, mC, muC, mu, Q] computes the massless to",
	"tal hadronic cross section.\"",
	(const char*)0,
	"RhadMass::usage = \"RhadMass[scheme, current, orderAlpha, runAlph",
	"a, runMass, order, nf, Mz, GammaZ, sin2ThetaW, aMz, mT, muT, mB,",
	" muB, mC, muC, mu, Q] computes the massless total hadronic cross",
	" section.\"",
	(const char*)0,
	"LambdaQCD::usage = \"LambdaQCD[scheme, order, runAlpha, run, nf, ",
	"Mz, aMz, mT, muT, mB, muB, mC, muC, mu] computes the running of ",
	"the quark masses with flavor matching.\"",
	(const char*)0,
	"Hyper2F1::usage=\"Hyper2F1[a, b, c, x] Hypergeometric Function in",
	" Fortran\"",
	(const char*)0,
	"HyperF32Exact::usage=\"HyperF32Exact[w, x] Hypergeometric Functio",
	"n in Fortran\"",
	(const char*)0,
	"InteCorre::usage=\"InteCorre[b, x0, x1] Incomplete Gamma Function",
	"\"",
	(const char*)0,
	"DiffDeltaGap::usage = \"DiffDeltaGap[gap, scheme, order, R0, R1, ",
	"mu0, mu1, muLambda, orderAlpha, runAlpha, nf, Mz, aMz, mT, muT, ",
	"mB, muB, mC, muC] computes the running of the gap parameter.\"",
	(const char*)0,
	"DiffDeltaGapMass::usage = \"DiffDeltaGapMass[gap, order, R0, R1, ",
	"mu0, mu1, muM, muLambda1, muLambda2, orderAlpha, runAlpha, runMa",
	"ss, nf, Mz, aMz, mT, muT, mB, muB, mC, muC] computes the running",
	" of the gap parameter with flavor matching.\"",
	(const char*)0,
	"Kernel::usage = \"Kernel[n, width, w, mu, p] computes the first n",
	"+1 kernels needed for resummation\"",
	(const char*)0,
	"GammaDerList::usage = \"GammaDerList[n, w] computes the first n+1",
	" derivatives of 1/Gamma\"",
	(const char*)0,
	"polyGamma::usage = \"polyGamma[n, w] computes the first n+1 deriv",
	"atives of Gamma\"",
	(const char*)0,
	"NGLKernel::usage = \"NGLKernel[n, n1, n2, width, w, mu, p] comput",
	"es the first 2*n kernels needed for NGL resummation\"",
	(const char*)0,
	"Taylor::usage = \"Taylor[c, lambda, k] computes the Taylor expans",
	"ion of the shape function\"",
	(const char*)0,
	"Model::usage = \"Model[c, lambda, k, l] computes the shape functi",
	"on\"",
	(const char*)0,
	"ModelUnstable::usage = \"ModelUnstable[shape, mt, Q, c, lambda, n",
	", k, l] computes the shape function convoluted with the unstable",
	" distribution\"",
	(const char*)0,
	"BreitModelUnstable::usage = \"BreitModelUnstable[shape, mt, Q, ga",
	"mma, c, lambda, n, k, l] computes the shape function convoluted ",
	"with the unstable distribution plus a BreitWigner\"",
	(const char*)0,
	"BreitModel::usage = \"BreitModel[c, lambda, width, k, l] computes",
	" the shape function convoluted with a Breit Wigner\"",
	(const char*)0,
	"MomentModel::usage = \"MomentModel[c, lambda, k] computes the sha",
	"pe function\"",
	(const char*)0,
	"ModelPiece::usage = \"ModelPiece[c, lambda, k, l] computes the sh",
	"ape function\"",
	(const char*)0,
	"TaylorPiece::usage = \"TaylorPiece[c, lambda, k] computes the Tay",
	"lor coefficients of the shape function\"",
	(const char*)0,
	"DeltaMSbar::usage = \"DeltaMSbar[order, runAlpha, run, nf, Mz, aM",
	"z, mT, muT, mB, muB, mC, muC, mu] computes the running of the qu",
	"ark masses with flavor matching.\"",
	(const char*)0,
	"JetMass::usage = \"JetMass[orderAlpha, runAlpha, order, run, nf, ",
	"Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, R, mu] computes th",
	"e Jet Mass running of the quark masses with flavor matching.\"",
	(const char*)0,
	"mmFromJetMass::usage = \"mmFromJetMass[orderAlpha, runAlpha, orde",
	"r, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, muLambda, R, mu]",
	" computes the Jet Mass running of the quark masses with flavor m",
	"atching.\"",
	(const char*)0,
	"Singular::usage = \"Singular[hard, shape, setup, gap, space, cum,",
	" orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, aMz, mT, ",
	"muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, lam",
	"bda, R0, mu0, delta0, h, tau] computes the Singular Thrust and C",
	"-parameter distrubution\"",
	(const char*)0,
	"SingularList::usage = \"SingularList[hard, shape, gap, space, cum",
	", orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, aMz, mT,",
	" muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, clen,",
	" lambda, R0, mu0, delta0, h, tau] computes the pieces of the Sin",
	"gular Thrust and C-parameter distrubution\"",
	(const char*)0,
	"MassNonDist::usage = \"MassNonDist[hard, shape, Eshape, setup, ga",
	"p, space, cum, scheme, orderAlpha, runAlpha, orderMass, runMass,",
	" order, run, nf, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambd",
	"a1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c, lambda, R",
	"0, mu0, delta0, h, tau] computes the Non distributional part of ",
	"the Singular Massive Thrust and C-parameter distrubution\"",
	(const char*)0,
	"MassNonDistPiece::usage = \"MassNonDistPiece[hard, shape, Eshape,",
	" gap, space, cum, scheme, orderAlpha, runAlpha, orderMass, runMa",
	"ss, order, run, nf, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLa",
	"mbda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu, c, lambda",
	", R0, mu0, delta0, h, tau] computes the Non distributional part ",
	"of the Singular Massive Thrust and C-parameter distrubution\"",
	(const char*)0,
	"SingularMass::usage = \"SingularMass[hard, shape, Eshape, setup, ",
	"gap, space, cum, scheme, abs, current, xi, xiB, orderAlpha, runA",
	"lpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ, amZ, m",
	"T, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS",
	", R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gamma",
	"Z, sin2ThetaW, tau] computes the Singular Massive Thrust and C-p",
	"arameter distrubution\"",
	(const char*)0,
	"SingularMassPiece::usage = \"SingularMassPiece[hard, shape, Eshap",
	"e, setup, gap, space, cum, scheme, abs, current, xi, xiB, orderA",
	"lpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, ",
	"mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH",
	", muJ, muS, R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0",
	", h, gammaZ, sin2ThetaW, tau] computes the Singular Massive Thru",
	"st and C-parameter distrubution\"",
	(const char*)0,
	"SingularPiece::usage = \"SingularPiece[hard, shape, gap, space, c",
	"um, orderAlpha, runAlpha, order, run, nf, j3, s3, G3, mZ, aMz, m",
	"T, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c, ",
	"lambda, R0, mu0, delta0, h, tau] computes the Singular Thrust an",
	"d C-parameter distrubution\"",
	(const char*)0,
	"SingularHJM::usage = \"SingularHJM[hard, setup, gap, space, cum, ",
	"orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32, G",
	"3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, mu",
	"S, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singu",
	"lar Thrust and C-parameter distrubution\"",
	(const char*)0,
	"SingularHJMPiece::usage = \"SingularHJMPiece[hard, gap, space, cu",
	"m, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, G3, mZ, ",
	"aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, m",
	"u, c, lambda, R0, mu0, delta0, h, tau] computes the Singular Thr",
	"ust and C-parameter distrubution\"",
	(const char*)0,
	"SingularDouble::usage = \"SingularDouble[hard, setup, gap, space,",
	" cum1, cum2, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3",
	", s31, s32, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, mu",
	"Lambda, Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, ",
	"rho1, rho2] computes the Singular Thrust and C-parameter distrub",
	"ution\"",
	(const char*)0,
	"SingularDoublePiece::usage = \"SingularDoublePiece[hard, gap, spa",
	"ce, cum1, cum2, orderAlpha, runAlpha, order, run, isoft, nf, j3,",
	" s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda,",
	" Q, muH, muJ, muS, R, mu, c, lambda, R0, mu0, delta0, h, rho1, r",
	"ho2] computes the Singular Thrust and C-parameter distrubution\"",
	(const char*)0,
	"SingularHJM1D::usage = \"SingularHJM1D[hard, gap, cum, orderAlpha",
	", runAlpha, order, run, isoft, nf, j3, s3, s31, s32, G3, mZ, aMz",
	", mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, ",
	"c, lambda, R0, mu0, delta0, h, tau] computes the Singular HJM di",
	"strubution with a 1D shape function\"",
	(const char*)0,
	"SingularHJM1DPiece::usage = \"SingularHJM1DPiece[hard, gap, cum, ",
	"orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32, G",
	"3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, mu",
	"S, R, mu, c, lambda, R0, mu0, delta0, h, tau] computes the Singu",
	"lar HJM distrubution with the piece of a 1D shape function\"",
	(const char*)0,
	"EWFactors::usage = \"EWFactors[nf, Q, Mz, GammaZ, sin2ThetaW] ele",
	"ctroweak factors\"",
	(const char*)0,
	"NGLIntegral::usage = \"NGLIntegral[nf, pow, w1, w2] computes the ",
	"exact NGL integral\"",
	(const char*)0,
	"NGLDoubleIntegral::usage = \"NGLDoubleIntegral[nf, pow, w1, w2, r",
	"] computes the exact NGL integral\"",
	(const char*)0,
	"ThrustNS1loop::usage = \"ThrustNS1loop[tau] computes the 1-loop t",
	"hrust NS function and its first two derivatives\"",
	(const char*)0,
	"FOMass::usage = \"FOMass[shape, current, m, Q, Mz, gammaZ, sin2Th",
	"etaW, tau] computes the 1-loop massive thrust NS function\"",
	(const char*)0,
	"ThrustNS2loop::usage = \"ThrustNS2loop[er, tau] computes the 2-lo",
	"op thrust NS function and its first two derivatives\"",
	(const char*)0,
	"NSMass::usage = \"NSMass[shape, setup, gap, cum, scheme, abs, cur",
	"rent, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, ",
	"Mz, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1, muLambda2, Q",
	", mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0, delta0",
	", h, gammaZ, sin2ThetaW, t] computes the NonSingular Massive Thr",
	"ust, HJM, SJM and C-parameter distribution, with a 1D model.\"",
	(const char*)0,
	"NSMassPiece::usage = \"NSMassPiece[shape, gap, cum, scheme, abs, ",
	"current, orderAlpha, runAlpha, orderMass, runMass, order, run, n",
	"f, Mz, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1, muLambda2",
	", Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0, del",
	"ta0, h, gammaZ, sin2ThetaW, t] computes the NonSingular Massive ",
	"Thrust, HJM, SJM, and C-parameter distribution, with a 1D model.",
	"\"",
	(const char*)0,
	"HJMNSMass::usage = \"HJMNSMass[setup, gap, cum, scheme, abs, curr",
	"ent, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, M",
	"z, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1, muLambda2, Q,",
	" mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0, delta0,",
	" h, gammaZ, sin2ThetaW, t] computes the NonSingular Massive HJM ",
	"distribution, with a 2-D model.\"",
	(const char*)0,
	"Begin[\"`Private`\"]",
	(const char*)0,
	"Print[\"You can access the complete function list typing '?Calipe",
	"r`*' \"]",
	(const char*)0,
	"realQ = Head[# + 1.] === Real &",
	(const char*)0,
	"End[]",
	(const char*)0,
	"EndPackage[]",
	(const char*)0,
	(const char*)0
};
#define CARDOF_EVALSTRS 111

static int _definepattern P(( MLINK, char*, char*, int));

static int _doevalstr P(( MLINK, int));

int  _MLDoCallPacket P(( MLINK, struct func[], int));


#if MLPROTOTYPES
int MLInstall( MLINK mlp)
#else
int MLInstall(mlp) MLINK mlp;
#endif
{
	int _res;
	_res = MLConnect(mlp);
	if (_res) _res = _doevalstr( mlp, 0);
	if (_res) _res = _doevalstr( mlp, 1);
	if (_res) _res = _doevalstr( mlp, 2);
	if (_res) _res = _doevalstr( mlp, 3);
	if (_res) _res = _doevalstr( mlp, 4);
	if (_res) _res = _doevalstr( mlp, 5);
	if (_res) _res = _doevalstr( mlp, 6);
	if (_res) _res = _doevalstr( mlp, 7);
	if (_res) _res = _doevalstr( mlp, 8);
	if (_res) _res = _doevalstr( mlp, 9);
	if (_res) _res = _doevalstr( mlp, 10);
	if (_res) _res = _doevalstr( mlp, 11);
	if (_res) _res = _doevalstr( mlp, 12);
	if (_res) _res = _doevalstr( mlp, 13);
	if (_res) _res = _doevalstr( mlp, 14);
	if (_res) _res = _doevalstr( mlp, 15);
	if (_res) _res = _doevalstr( mlp, 16);
	if (_res) _res = _doevalstr( mlp, 17);
	if (_res) _res = _doevalstr( mlp, 18);
	if (_res) _res = _doevalstr( mlp, 19);
	if (_res) _res = _doevalstr( mlp, 20);
	if (_res) _res = _doevalstr( mlp, 21);
	if (_res) _res = _doevalstr( mlp, 22);
	if (_res) _res = _doevalstr( mlp, 23);
	if (_res) _res = _doevalstr( mlp, 24);
	if (_res) _res = _doevalstr( mlp, 25);
	if (_res) _res = _doevalstr( mlp, 26);
	if (_res) _res = _doevalstr( mlp, 27);
	if (_res) _res = _doevalstr( mlp, 28);
	if (_res) _res = _doevalstr( mlp, 29);
	if (_res) _res = _doevalstr( mlp, 30);
	if (_res) _res = _doevalstr( mlp, 31);
	if (_res) _res = _doevalstr( mlp, 32);
	if (_res) _res = _doevalstr( mlp, 33);
	if (_res) _res = _doevalstr( mlp, 34);
	if (_res) _res = _doevalstr( mlp, 35);
	if (_res) _res = _doevalstr( mlp, 36);
	if (_res) _res = _doevalstr( mlp, 37);
	if (_res) _res = _doevalstr( mlp, 38);
	if (_res) _res = _doevalstr( mlp, 39);
	if (_res) _res = _doevalstr( mlp, 40);
	if (_res) _res = _doevalstr( mlp, 41);
	if (_res) _res = _doevalstr( mlp, 42);
	if (_res) _res = _doevalstr( mlp, 43);
	if (_res) _res = _doevalstr( mlp, 44);
	if (_res) _res = _doevalstr( mlp, 45);
	if (_res) _res = _doevalstr( mlp, 46);
	if (_res) _res = _doevalstr( mlp, 47);
	if (_res) _res = _doevalstr( mlp, 48);
	if (_res) _res = _doevalstr( mlp, 49);
	if (_res) _res = _doevalstr( mlp, 50);
	if (_res) _res = _doevalstr( mlp, 51);
	if (_res) _res = _doevalstr( mlp, 52);
	if (_res) _res = _doevalstr( mlp, 53);
	if (_res) _res = _doevalstr( mlp, 54);
	if (_res) _res = _doevalstr( mlp, 55);
	if (_res) _res = _doevalstr( mlp, 56);
	if (_res) _res = _doevalstr( mlp, 57);
	if (_res) _res = _doevalstr( mlp, 58);
	if (_res) _res = _doevalstr( mlp, 59);
	if (_res) _res = _doevalstr( mlp, 60);
	if (_res) _res = _doevalstr( mlp, 61);
	if (_res) _res = _doevalstr( mlp, 62);
	if (_res) _res = _doevalstr( mlp, 63);
	if (_res) _res = _doevalstr( mlp, 64);
	if (_res) _res = _doevalstr( mlp, 65);
	if (_res) _res = _doevalstr( mlp, 66);
	if (_res) _res = _doevalstr( mlp, 67);
	if (_res) _res = _doevalstr( mlp, 68);
	if (_res) _res = _doevalstr( mlp, 69);
	if (_res) _res = _doevalstr( mlp, 70);
	if (_res) _res = _doevalstr( mlp, 71);
	if (_res) _res = _doevalstr( mlp, 72);
	if (_res) _res = _doevalstr( mlp, 73);
	if (_res) _res = _doevalstr( mlp, 74);
	if (_res) _res = _doevalstr( mlp, 75);
	if (_res) _res = _doevalstr( mlp, 76);
	if (_res) _res = _doevalstr( mlp, 77);
	if (_res) _res = _doevalstr( mlp, 78);
	if (_res) _res = _doevalstr( mlp, 79);
	if (_res) _res = _doevalstr( mlp, 80);
	if (_res) _res = _doevalstr( mlp, 81);
	if (_res) _res = _doevalstr( mlp, 82);
	if (_res) _res = _doevalstr( mlp, 83);
	if (_res) _res = _doevalstr( mlp, 84);
	if (_res) _res = _doevalstr( mlp, 85);
	if (_res) _res = _doevalstr( mlp, 86);
	if (_res) _res = _doevalstr( mlp, 87);
	if (_res) _res = _doevalstr( mlp, 88);
	if (_res) _res = _doevalstr( mlp, 89);
	if (_res) _res = _doevalstr( mlp, 90);
	if (_res) _res = _doevalstr( mlp, 91);
	if (_res) _res = _doevalstr( mlp, 92);
	if (_res) _res = _doevalstr( mlp, 93);
	if (_res) _res = _doevalstr( mlp, 94);
	if (_res) _res = _doevalstr( mlp, 95);
	if (_res) _res = _doevalstr( mlp, 96);
	if (_res) _res = _doevalstr( mlp, 97);
	if (_res) _res = _doevalstr( mlp, 98);
	if (_res) _res = _doevalstr( mlp, 99);
	if (_res) _res = _doevalstr( mlp, 100);
	if (_res) _res = _doevalstr( mlp, 101);
	if (_res) _res = _doevalstr( mlp, 102);
	if (_res) _res = _doevalstr( mlp, 103);
	if (_res) _res = _doevalstr( mlp, 104);
	if (_res) _res = _doevalstr( mlp, 105);
	if (_res) _res = _doevalstr( mlp, 106);
	if (_res) _res = _doevalstr( mlp, 107);
	if (_res) _res = _definepattern(mlp, (char *)"EWFactors[nf_, Q_, Mz_, GammaZ_, sin2ThetaW_]", (char *)"{nf, Q, Mz, GammaZ, sin2ThetaW}", 0);
	if (_res) _res = _definepattern(mlp, (char *)"LegendreList[n_, k_, x_]", (char *)"{n, k, x}", 1);
	if (_res) _res = _definepattern(mlp, (char *)"QLegendreList[n_, x_]", (char *)"{n, x}", 2);
	if (_res) _res = _definepattern(mlp, (char *)"GammaR[str_, nf_]", (char *)"{str, nf}", 3);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessProf[terms_, hard_, shape_, setup_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,                 delta0_, h_, tau_]", (char *)"{terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,                  R0, muR0, delta0, h, tau}", 4);
	if (_res) _res = _definepattern(mlp, (char *)"FindOrigin[shape_, gap_, orderAlpha_, runAlpha_, order_, run_, nf_,                 mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_,                 n0_, n1_, t2_, tR_, ts_, slope_, cnt_, eH_, eS_, eR_, R0_, muR0_, delta0_,                 h_]", (char *)"{shape, gap, orderAlpha, runAlpha, order, run, nf, mZ, amZ, mT, muT, mB,                  muB, mC, muC, muLambda, Q, mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH,                  eS, eR, R0, muR0, delta0, h}", 5);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessProfPiece[terms_, hard_, shape_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,                 muR0_, delta0_, h_, tau_]", (char *)"{terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,                  lambda, R0, muR0, delta0, h, tau}", 6);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessProfPieceList[terms_, hard_, shape_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,                 muR0_, delta0_, h_, tauList_]", (char *)"{terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,                  lambda, R0, muR0, delta0, h, tauList}", 7);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessPieceBin[terms_, hard_, shape_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,                 muR0_, delta0_, h_, tauList_]", (char *)"{terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,                  lambda, R0, muR0, delta0, h, Flatten[tauList], Length[tauList]}", 8);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessProfPiece[terms_, hard_, shape_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, clen_, lambda_, R0_,                 muR0_, delta0_, h_, tau_, tau2_]", (char *)"{terms, hard, shape, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, clen,                  lambda, R0, muR0, delta0, h, tau, tau2}", 9);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessProfList[terms_, hard_, shape_, setup_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,                 delta0_, h_, tauList_]", (char *)"{terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,                  R0, muR0, delta0, h, tauList}", 10);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessBinList[terms_, hard_, shape_, setup_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,                 delta0_, h_, tauList_]", (char *)"{terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,                  R0, muR0, delta0, h, Flatten[tauList], Length[tauList]}", 11);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessProf[terms_, hard_, shape_, setup_, gap_, space_, cum_,                 orderAlpha_, runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_,                 tR_, ts_, slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_,                 delta0_, h_, tau_, tau2_]", (char *)"{terms, hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order,                  run, nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q,                  mu0, Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda,                  R0, muR0, delta0, h, tau, tau2}", 12);
	if (_res) _res = _definepattern(mlp, (char *)"MasslessMoment[terms_, hard_, shape_, setup_, gap_, space_, orderAlpha_,                 runAlpha_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda_, Q_, mu0_, Rat0_, n0_, n1_, t2_, tR_, ts_,                 slope_, cnt_, eH_, eS_, eJ_, eR_, ns_, c_, lambda_, R0_, muR0_, delta0_,                 h_, tau_, tau2_, pow_]", (char *)"{terms, hard, shape, setup, gap, space, orderAlpha, runAlpha, order, run,                  nf, j3, s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda, Q, mu0,                  Rat0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, c, lambda, R0,                  muR0, delta0, h, tau, tau2, pow}", 13);
	if (_res) _res = _definepattern(mlp, (char *)"Profiles[Q_, mu0_, R0_, n0_, n1_, t2_, tR_, ts_, slope_, cnt_, eH_, eS_,                 eJ_, eR_, ns_, tau_]", (char *)"{Q, mu0, R0, n0, n1, t2, tR, ts, slope, cnt, eH, eS, eJ, eR, ns, tau}", 14);
	if (_res) _res = _definepattern(mlp, (char *)"ProfilesMass[Q_, beta_, mu0_, delLamb_, R0_, n0_, delta0_, n1_,                 delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_, mass_, muM_, ns_,                 def_, EShape_, tau_]", (char *)"{Q, beta, mu0, delLamb, R0, n0, delta0, n1, delta1, t2, ts, slope,                  cnt, eH, eS, eJ, mass, muM, ns, def, EShape, tau}", 15);
	if (_res) _res = _definepattern(mlp, (char *)"MCtop[shape_, mt_, Q_, n_, k_, x_]", (char *)"{shape, mt, Q, n, k, x}", 16);
	if (_res) _res = _definepattern(mlp, (char *)"BreitUnstable[shape_, mt_, Q_, gamma_, n_, k_, x_]", (char *)"{shape, mt, Q, gamma, n, k, x}", 17);
	if (_res) _res = _definepattern(mlp, (char *)"DeltaMCtop[shape_, mt_, Q_]", (char *)"{shape, mt, Q}", 18);
	if (_res) _res = _definepattern(mlp, (char *)"MassIntersection[Q_, beta_, mu0_, delLamb_, n0_, delta0_, n1_,delta1_,                 t2_, ts_, slope_, cnt_, eH_, eS_, eJ_, mass_, muM_, def_, EShape_]", (char *)"{Q, beta, mu0, delLamb, n0, delta0, n1, delta1, t2, ts, slope,                  cnt, eH, eS, eJ, mass, muM, def, EShape}", 19);
	if (_res) _res = _definepattern(mlp, (char *)"NSMassPiece[shape_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,                 runAlpha_, order_, run_, orderMass_ runMass_, nf_, mZ_, aMz_, mT_, muT_,                 mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,                 muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,                 sin2ThetaW_, tau_]", (char *)"{shape, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,                 orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,                 muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0,                 delta0, h, gammaZ, sin2ThetaW, tau}", 20);
	if (_res) _res = _definepattern(mlp, (char *)"NSMassPiece[shape_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,                 runAlpha_, order_, run_, orderMass_ runMass_, nf_, mZ_, aMz_, mT_, muT_,                 mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,                 muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,                 sin2ThetaW_, tau_, tau2_]", (char *)"{shape, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,                 orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,                 muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c, lambda, R0, mu0,                 delta0, h, gammaZ, sin2ThetaW, tau, tau2}", 21);
	if (_res) _res = _definepattern(mlp, (char *)"NSMass[shape_, setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,                 runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,                 mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,                 muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,                 sin2ThetaW_, tau_]", (char *)"{shape, setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha,                 order, run, orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC,                 muC, muLambda1, muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c,                 lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau}", 22);
	if (_res) _res = _definepattern(mlp, (char *)"NSMass[shape_, setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,                 runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,                 mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,                 muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,                 sin2ThetaW_, tau_, tau2_]", (char *)"{shape, setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha,                 order, run, orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC,                 muC, muLambda1, muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, c,                 lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2}", 23);
	if (_res) _res = _definepattern(mlp, (char *)"HJMNSMass[setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,                 runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,                 mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,                 muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,                 sin2ThetaW_, tau_]", (char *)"{setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,                 orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,                 muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, Flatten[Transpose[c]],                 Length[c], lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau}", 24);
	if (_res) _res = _definepattern(mlp, (char *)"HJMNSMass[setup_, gap_, cum_, scheme_, abs_, current_, orderAlpha_,                 runAlpha_, order_, run_, orderMass_, runMass_, nf_, mZ_, aMz_, mT_, muT_,                 mB_, muBottom_, mC_, muC_, muLambda1_, muLambda2_, Q_, mu_, muM_, muB_,                 muS_, R_, Rmass_, width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_,                 sin2ThetaW_, tau_, tau2_]", (char *)"{setup, gap, cum, scheme, abs, current, orderAlpha, runAlpha, order, run,                 orderMass, runMass, nf, mZ, aMz, mT, muT, mB, muBottom, mC, muC, muLambda1,                 muLambda2, Q, mu, muM, muB, muS, R, Rmass, width, Flatten[Transpose[c]],                 Length[c], lambda, R0, mu0, delta0, h, gammaZ, sin2ThetaW, tau, tau2}", 25);
	if (_res) _res = _definepattern(mlp, (char *)"Singular[hard_, shape_, setup_, gap_, space_, cum_, orderAlpha_, runAlpha_,                 order_, run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_,                 muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_,                 delta0_, h_, tau_]", (char *)"{hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3,                  s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,                  R, mu, c, lambda, R0, mu0, delta0, h, tau}", 26);
	if (_res) _res = _definepattern(mlp, (char *)"SingularList[hard_, shape_, gap_, space_, cum_, orderAlpha_, runAlpha_,                 order_, run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_,                 muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, clen_, lambda_, R0_, mu0_,                 delta0_, h_, tau_]", (char *)"{hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3,                  s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,                  R, mu, clen, lambda, R0, mu0, delta0, h, tau}", 27);
	if (_res) _res = _definepattern(mlp, (char *)"Singular[hard_, shape_, setup_, gap_, space_, cum_, orderAlpha_, runAlpha_,                 order_, run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_,                 muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_,                 delta0_, h_, tau1_, tau2_]", (char *)"{hard, shape, setup, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3,                  s3, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,                  R, mu, c, lambda, R0, mu0, delta0, h, tau1, tau2}", 28);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveProf[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tau_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, tau}", 29);
	if (_res) _res = _definepattern(mlp, (char *)"MassOrigin[shape_, Eshape_, gap_, scheme_, orderAlpha_, runAlpha_,                 orderMass_, runMass_, order_, run_, nf_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, R0_, muR0_, del0_, h_]", (char *)"{shape, Eshape, gap, scheme, orderAlpha, runAlpha, orderMass, runMass,                  order, run, nf, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2,                  Q, beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope,                  cnt, eH, eS, eJ, mass, muM, R0, muR0, del0, h}", 30);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveProfPiece[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tau_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, tau}", 31);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveProfPieceList[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tauList_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, tauList}", 32);
	if (_res) _res = _definepattern(mlp, (char *)"MassivePieceBin[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tauList_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, Flatten[tauList], Length[tauList]}", 33);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveProfPiece[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, clen_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tau_, tau2_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, clen, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, tau, tau2}", 34);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveProf[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tau_, tau2_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, tau, tau2}", 35);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveProfList[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tauList_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, tauList}", 36);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveBinList[terms_, hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, del0_, h_, gammaZ_,                 sin2ThetaW_, tauList_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, muR0, del0, h, gammaZ,                  sin2ThetaW, Flatten[tauList], Length[tauList]}", 37);
	if (_res) _res = _definepattern(mlp, (char *)"MassiveMoment[terms_, hard_, shape_, Eshape_, setup_, gap_, space_,                 scheme_, abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_,                 runMass_, order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda1_, muLambda2_, Q_, beta_, mu0_, deltaLambda_,                 Rat0_, n0_, delta0_, n1_, delta1_, t2_, ts_, slope_, cnt_, eH_, eS_, eJ_,                 mass_, muM_, ns_, width_, c_, lambda_, R0_, muR0_, delta0_, h_, gammaZ_,                 sin2ThetaW_, tau_, tau2_, pow_]", (char *)"{terms, hard, shape, Eshape, setup, gap, space, scheme, abs, current,                  xi, xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3,                  s3, G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q,                  beta, mu0, deltaLambda, Rat0, n0, delta0, n1, delta1, t2, ts, slope, cnt,                  eH, eS, eJ, mass, muM, ns, width, c, lambda, R0, mu0, delta0, h, gammaZ,                  sin2ThetaW, tau, tau2, pow}", 38);
	if (_res) _res = _definepattern(mlp, (char *)"SingularMass[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,                 abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,                 order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,                 muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,                 width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_]", (char *)"{hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi,                  xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3,                  mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ,                  muS, R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ,                  sin2ThetaW, tau}", 39);
	if (_res) _res = _definepattern(mlp, (char *)"SingularMass[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,                 abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,                 order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,                 muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,                 width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_,                 tau2_]", (char *)"{hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi,                  xiB, orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3,                  G3, mZ, amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH,                  muJ, muS, R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h,                  gammaZ, sin2ThetaW, tau, tau2}", 40);
	if (_res) _res = _definepattern(mlp, (char *)"MassNonDist[hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, orderAlpha_, runAlpha_, orderMass_, runMass_, order_,                 run_, nf_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_,                 muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_,                 mu_, c_, lambda_, R0_, mu0_, delta0_, h_, tau_]", (char *)"{hard, shape, Eshape, setup, gap, space, cum, scheme, orderAlpha,                  runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,                  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,                  muM, mu, c, lambda, R0, mu0, delta0, h, tau}", 41);
	if (_res) _res = _definepattern(mlp, (char *)"MassNonDist[hard_, shape_, Eshape_, setup_, gap_, space_, cum_,                 scheme_, orderAlpha_, runAlpha_, orderMass_, runMass_, order_,                 run_, nf_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_,                 muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_,                 mu_, c_, lambda_, R0_, mu0_, delta0_, h_, tau_, tau2_]", (char *)"{hard, shape, Eshape, setup, gap, space, cum, scheme, orderAlpha,                  runAlpha, orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT,                  mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass,                  muM, mu, c, lambda, R0, mu0, delta0, h, tau, tau2}", 42);
	if (_res) _res = _definepattern(mlp, (char *)"MassNonDistPiece[hard_, shape_, Eshape_, gap_, space_, cum_, scheme_,                  orderAlpha_, runAlpha_, orderMass_, runMass_, order_, run_, nf_, G3_,                  mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_, muLambda1_, muLambda2_,                  Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_, c_, lambda_, R0_, mu0_,                  delta0_, h_, tau_]", (char *)"{hard, shape, Eshape, gap, space, cum, scheme, orderAlpha, runAlpha,                  orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, mB, muB,                  mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,                  c, lambda, R0, mu0, delta0, h, tau}", 43);
	if (_res) _res = _definepattern(mlp, (char *)"MassNonDistPiece[hard_, shape_, Eshape_, gap_, space_, cum_, scheme_,                  orderAlpha_, runAlpha_, orderMass_, runMass_, order_, run_, nf_, G3_,                  mZ_, amZ_, mT_, muT_, mB_, muB_, mC_, muC_, muLambda1_, muLambda2_,                  Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_, c_, lambda_, R0_, mu0_,                  delta0_, h_, tau_, tau2_]", (char *)"{hard, shape, Eshape, gap, space, cum, scheme, orderAlpha, runAlpha,                  orderMass, runMass, order, run, nf, G3, mZ, amZ, mT, muT, mB, muB,                  mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS, R, Rmass, muM, mu,                  c, lambda, R0, mu0, delta0, h, tau, tau2}", 44);
	if (_res) _res = _definepattern(mlp, (char *)"SingularMassPiece[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,                 abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,                 order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,                 muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,                 width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_]", (char *)"{hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi, xiB,                  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ,                  amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS,                  R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ,                  sin2ThetaW, tau}", 45);
	if (_res) _res = _definepattern(mlp, (char *)"SingularMassPiece[hard_, shape_, Eshape_, setup_, gap_, space_, cum_, scheme_,                 abs_, current_, xi_, xiB_, orderAlpha_, runAlpha_, orderMass_, runMass_,                 order_, run_, nf_, j3_, s3_, G3_, mZ_, amZ_, mT_, muT_, mB_, muB_, mC_,                 muC_, muLambda1_, muLambda2_, Q_, muH_, muJ_, muS_, R_, Rmass_, muM_, mu_,                 width_, c_, lambda_, R0_, mu0_, delta0_, h_, gammaZ_, sin2ThetaW_, tau_,                 tau2_]", (char *)"{hard, shape, Eshape, setup, gap, space, cum, scheme, abs, current, xi, xiB,                  orderAlpha, runAlpha, orderMass, runMass, order, run, nf, j3, s3, G3, mZ,                  amZ, mT, muT, mB, muB, mC, muC, muLambda1, muLambda2, Q, muH, muJ, muS,                  R, Rmass, muM, mu, width, c, lambda, R0, mu0, delta0, h, gammaZ,                  sin2ThetaW, tau, tau2}", 46);
	if (_res) _res = _definepattern(mlp, (char *)"SingularHJM[hard_, setup_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_,                 run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,                 R0_, mu0_, delta0_, h_, tau_]", (char *)"{hard, setup, gap, space, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3,                  s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH,                  muJ, muS, R, mu, Flatten[Transpose[c]], Length[c], lambda, R0, mu0,                  delta0, h, tau}", 47);
	if (_res) _res = _definepattern(mlp, (char *)"SingularHJM1D[hard_, gap_, cum_, orderAlpha_, runAlpha_, order_,                 run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_,                 muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,                 R0_, mu0_, delta0_, h_, tau_]", (char *)"{hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31, s32,                  G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,                  R, mu, c, lambda, R0, mu0, delta0, h, tau}", 48);
	if (_res) _res = _definepattern(mlp, (char *)"SingularHJM1DPiece[hard_, gap_, cum_, orderAlpha_, runAlpha_, order_, run_,                 isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_,                 mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_,                 mu0_, delta0_, h_, tau_]", (char *)"{hard, gap, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3,                  s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,                  R, mu, c, lambda, R0, mu0, delta0, h, tau}", 49);
	if (_res) _res = _definepattern(mlp, (char *)"SingularDouble[hard_, setup_, gap_, space_, cum1_, cum2_, orderAlpha_, runAlpha_,                 order_, run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,                 R0_, mu0_, delta0_, h_, rho1_, rho2_]", (char *)"{hard, setup, gap, space, cum1, cum2, orderAlpha, runAlpha, order, run, isoft,                  nf, j3, s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH,                  muJ, muS, R, mu, Flatten[Transpose[c]], Length[c], lambda, R0, mu0,                  delta0, h, rho1, rho2}", 50);
	if (_res) _res = _definepattern(mlp, (char *)"SingularHJMPiece[hard_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_, run_,                 isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_,                 muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_, delta0_,                 h_, tau_]", (char *)"{hard, gap, space, cum, orderAlpha, runAlpha, order, run, isoft, nf, j3, s3, s31,                  s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R,                  mu, c, lambda, R0, mu0, delta0, h, tau}", 51);
	if (_res) _res = _definepattern(mlp, (char *)"SingularDoublePiece[hard_, gap_, space_, cum1_, cum2_, orderAlpha_, runAlpha_,                 order_, run_, isoft_, nf_, j3_, s3_, s31_, s32_, G3_, mZ_, aMz_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_,                 R0_, mu0_, delta0_, h_, rho1_, rho2_]", (char *)"{hard, gap, space, cum1, cum2, orderAlpha, runAlpha, order, run, isoft, nf, j3,                  s3, s31, s32, G3, mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS,                  R, mu, c, lambda, R0, mu0, delta0, h, rho1, rho2}", 52);
	if (_res) _res = _definepattern(mlp, (char *)"SingularPiece[hard_, shape_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_,                 run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_,                 muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_, delta0_,                 h_, tau_]", (char *)"{hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3,                  mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c,                  lambda, R0, mu0, delta0, h, tau}", 53);
	if (_res) _res = _definepattern(mlp, (char *)"SingularPiece[hard_, shape_, gap_, space_, cum_, orderAlpha_, runAlpha_, order_,                 run_, nf_, j3_, s3_, G3_, mZ_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_,                 muLambda_, Q_, muH_, muJ_, muS_, R_, mu_, c_, lambda_, R0_, mu0_, delta0_,                 h_, tau_, tau2_]", (char *)"{hard, shape, gap, space, cum, orderAlpha, runAlpha, order, run, nf, j3, s3, G3,                  mZ, aMz, mT, muT, mB, muB, mC, muC, muLambda, Q, muH, muJ, muS, R, mu, c,                  lambda, R0, mu0, delta0, h, tau, tau2}", 54);
	if (_res) _res = _definepattern(mlp, (char *)"DiffDeltaGap[str_, scheme_, order_, R0_, R1_, mu0_, mu1_, muLambda_,                 orderAlpha_, runAlpha_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_]", (char *)"{str, scheme, order, R0, R1, mu0, mu1, muLambda, orderAlpha, runAlpha,                  nf, Mz, aMz, mT, muT, mB, muB, mC, muC}", 55);
	if (_res) _res = _definepattern(mlp, (char *)"DiffDeltaGapMass[str_, order_, R0_, R1_, mu0_, mu1_, muM_, muLambda1_,                 muLambda2_, orderAlpha_, runAlpha_, runMass_, nf_, Mz_, aMz_, mT_, muT_,                 mB_, muB_, mC_, muC_]", (char *)"{str, order, R0, R1, mu0, mu1, muM, muLambda1, muLambda2, orderAlpha,                  runAlpha, runMass, nf, Mz, aMz, mT, muT, mB, muB, mC, muC}", 56);
	if (_res) _res = _definepattern(mlp, (char *)"HyperF32Exact[w_, x_]", (char *)"{w, x}", 57);
	if (_res) _res = _definepattern(mlp, (char *)"Model[c_, lambda_, k_, l_]", (char *)"{c, lambda, k, l}", 58);
	if (_res) _res = _definepattern(mlp, (char *)"ModelUnstable[shape_, mt_, Q_, c_, lambda_, n_, k_, l_]", (char *)"{shape, mt, Q, c, lambda, n, k, l}", 59);
	if (_res) _res = _definepattern(mlp, (char *)"BreitModelUnstable[shape_, mt_, Q_, gamma_, c_, lambda_, n_, k_, l_]", (char *)"{shape, mt, Q, gamma, c, lambda, n, k, l}", 60);
	if (_res) _res = _definepattern(mlp, (char *)"ModelUnstable[shape_, mt_, Q_, c_, lambda_, n_, k_, l_, l2_]", (char *)"{shape, mt, Q, c, lambda, n, k, l, l2}", 61);
	if (_res) _res = _definepattern(mlp, (char *)"Model[c_, lambda_, k_, l_, l2_]", (char *)"{c, lambda, k, l, l2}", 62);
	if (_res) _res = _definepattern(mlp, (char *)"BreitModel[c_, lambda_, width_, k_, l_]", (char *)"{c, lambda, width, k, l}", 63);
	if (_res) _res = _definepattern(mlp, (char *)"BreitModel[c_, lambda_, width_, k_, l_, l2_]", (char *)"{c, lambda, width, k, l, l2}", 64);
	if (_res) _res = _definepattern(mlp, (char *)"Taylor[c_, lambda_, k_]", (char *)"{c, lambda, k}", 65);
	if (_res) _res = _definepattern(mlp, (char *)"MomentModel[c_, lambda_, k_]", (char *)"{c, lambda, k}", 66);
	if (_res) _res = _definepattern(mlp, (char *)"ModelPiece[c_, lambda_, k_, l_]", (char *)"{c, lambda, k, l}", 67);
	if (_res) _res = _definepattern(mlp, (char *)"TaylorPiece[c_, lambda_, k_]", (char *)"{c, lambda, k}", 68);
	if (_res) _res = _definepattern(mlp, (char *)"Hyper2F1[a_, b_, c_, x_]", (char *)"{a, b, c, x}", 69);
	if (_res) _res = _definepattern(mlp, (char *)"InteCorre[b_, x0_, x1_]", (char *)"{b, x0, x1}", 70);
	if (_res) _res = _definepattern(mlp, (char *)"AnomDim[str_, nf_, G4_]", (char *)"{str, nf, G4}", 71);
	if (_res) _res = _definepattern(mlp, (char *)"Scoef[str_, nf_]", (char *)"{str, nf}", 72);
	if (_res) _res = _definepattern(mlp, (char *)"Polylog[n_, z_]", (char *)"{n, z}", 73);
	if (_res) _res = _definepattern(mlp, (char *)"DiLog[z_]", (char *)"{z}", 74);
	if (_res) _res = _definepattern(mlp, (char *)"pFq[a_, b_, z_]", (char *)"{a, b, z}", 75);
	if (_res) _res = _definepattern(mlp, (char *)"Elliptic3[psi_, k_, c_]", (char *)"{psi, k, c}", 76);
	if (_res) _res = _definepattern(mlp, (char *)"NGLFunction[n_, z_]", (char *)"{n, z}", 77);
	if (_res) _res = _definepattern(mlp, (char *)"NGLSoft[n_, z_]", (char *)"{n, Re[z], Im[z]}", 78);
	if (_res) _res = _definepattern(mlp, (char *)"ComplexPolylog[n_, z_]", (char *)"{n, Re[z], Im[z]}", 79);
	if (_res) _res = _definepattern(mlp, (char *)"CLi2[z_]", (char *)"{Re[z], Im[z]}", 80);
	if (_res) _res = _definepattern(mlp, (char *)"ThrustNS1loop[t_]", (char *)"{t}", 81);
	if (_res) _res = _definepattern(mlp, (char *)"FOMass[shape_, current_, m_, Q_, Mz_, gammaZ_, sin2ThetaW_, t_]", (char *)"{shape, current, m, Q, Mz, gammaZ, sin2ThetaW, t}", 82);
	if (_res) _res = _definepattern(mlp, (char *)"ThrustNS2loop[er_, t_]", (char *)"{er, t}", 83);
	if (_res) _res = _definepattern(mlp, (char *)"CLi3[z_]", (char *)"{Re[z], Im[z]}", 84);
	if (_res) _res = _definepattern(mlp, (char *)"Delta[str_, nf_, mu_, R_]", (char *)"{str, nf, mu, R}", 85);
	if (_res) _res = _definepattern(mlp, (char *)"DeltaGap[str_, orderAlpha_, runAlpha_, runMass_, nf_, mZ_, aMz_, mT_,                 muT_, mB_, muB_, mC_, muC_, mu_, R_]", (char *)"{str, orderAlpha, runAlpha, runMass, nf, mZ, aMz, mT, muT, mB, muB, mC,                  muC, mu, R}", 86);
	if (_res) _res = _definepattern(mlp, (char *)"CoefMat[str_, nf_, s3_]", (char *)"{str, nf, s3}", 87);
	if (_res) _res = _definepattern(mlp, (char *)"wTilde[order_, nf_, gamma_, a0_, a1_]", (char *)"{order, nf, gamma, a0, a1}", 88);
	if (_res) _res = _definepattern(mlp, (char *)"kTilde[order_, nf_, gamma_, a0_, a1_]", (char *)"{order, nf, gamma, a0, a1}", 89);
	if (_res) _res = _definepattern(mlp, (char *)"AlphaQCD[str_, method_, order_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_, mC_,                  muC_, mu_]", (char *)"{str, method, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}", 90);
	if (_res) _res = _definepattern(mlp, (char *)"AlphaComplex[str_, method_, order_, run_, nf_, Mz_, aMz_, mT_,                 muT_, mB_, muB_, mC_, muC_, mu_]", (char *)"{str, method, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC,                 muC, Re[mu], Im[mu]}", 91);
	if (_res) _res = _definepattern(mlp, (char *)"MSbarMass[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_,                 mC_, muC_, mu_]", (char *)"{order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}", 92);
	if (_res) _res = _definepattern(mlp, (char *)"PoleMass[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_, mT_, muT_,                 mB_, muB_, mC_, muC_, mu_]", (char *)"{orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC,                  mu}", 93);
	if (_res) _res = _definepattern(mlp, (char *)"MSbarMassLow[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_,                 muB_, mC_, muC_, mu_]", (char *)"{order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}", 94);
	if (_res) _res = _definepattern(mlp, (char *)"MSRMass[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_,                 mC_, muC_, mu_, R_]", (char *)"{order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu, R}", 95);
	if (_res) _res = _definepattern(mlp, (char *)"mmfromMSR[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_,                 mC_, muC_, mu_, R_]", (char *)"{order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu, R}", 96);
	if (_res) _res = _definepattern(mlp, (char *)"mmfromMSRNatural[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_,                 mT_, muT_, mB_, muB_, mC_, muC_, mu_, R_]", (char *)"{orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC,                  mu, R}", 97);
	if (_res) _res = _definepattern(mlp, (char *)"MSRNaturalMass[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_,                 mT_, muT_, mB_, muB_, mC_, muC_, mu_, R_]", (char *)"{orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC,                  mu, R}", 98);
	if (_res) _res = _definepattern(mlp, (char *)"JetMass[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, R_, mu_]", (char *)"{orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC,                  muC, muLambda, R, mu}", 99);
	if (_res) _res = _definepattern(mlp, (char *)"mmFromJetMass[orderAlpha_, runAlpha_, order_, run_, nf_, Mz_, aMz_, mT_,                 muT_, mB_, muB_, mC_, muC_, muLambda_, R_, mu_]", (char *)"{orderAlpha, runAlpha, order, run, nf, Mz, aMz, mT, muT, mB, muB, mC,                  muC, muLambda, R, mu}", 100);
	if (_res) _res = _definepattern(mlp, (char *)"DeltaMSbar[order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_, mB_, muB_,                  mC_, muC_, mu_]", (char *)"{order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}", 101);
	if (_res) _res = _definepattern(mlp, (char *)"Rhad[scheme_, orderAlpha_, runAlpha_, order_, nf_, Mz_, aMz_, mT_, muT_,                  mB_, muB_, mC_, muC_, mu_, Q_]", (char *)"{scheme, orderAlpha, runAlpha, order, nf, Mz, aMz, mT, muT, mB, muB, mC,                  muC, mu, Q}", 102);
	if (_res) _res = _definepattern(mlp, (char *)"RhadMass[scheme_, current_, orderAlpha_, runAlpha_, runMass_, order_, nf_,                 Mz_, GammaZ_, sin2ThetaW_, aMz_, mT_, muT_, mB_, muB_, mC_, muC_, mu_, Q_]", (char *)"{scheme, current, orderAlpha, runAlpha, runMass, order, nf, Mz, GammaZ,                 sin2ThetaW, aMz, mT, muT, mB, muB, mC, muC, mu, Q}", 103);
	if (_res) _res = _definepattern(mlp, (char *)"LambdaQCD[scheme_, order_, runAlpha_, run_, nf_, Mz_, aMz_, mT_, muT_,                  mB_, muB_, mC_, muC_, mu_]", (char *)"{scheme, order, runAlpha, run, nf, Mz, aMz, mT, muT, mB, muB, mC, muC, mu}", 104);
	if (_res) _res = _definepattern(mlp, (char *)"Kernel[n_, width_, w_, mu_, p_]", (char *)"{n, width, w, mu, p}", 105);
	if (_res) _res = _definepattern(mlp, (char *)"GammaDerList[n_, w_]", (char *)"{n, w}", 106);
	if (_res) _res = _definepattern(mlp, (char *)"polyGamma[n_, w_]", (char *)"{n, w}", 107);
	if (_res) _res = _definepattern(mlp, (char *)"NGLKernel[n_, n1_, n2_, width_, w_, mu_, p_]", (char *)"{n, n1, n2, width, w, mu, p}", 108);
	if (_res) _res = _definepattern(mlp, (char *)"NGLIntegral[nf_, pow_, w1_, w2_]", (char *)"{nf, pow, w1, w2}", 109);
	if (_res) _res = _definepattern(mlp, (char *)"NGLDoubleIntegral[nf_, pow_, w1_, w2_, r_]", (char *)"{nf, pow, w1, w2, r}", 110);
	if (_res) _res = _doevalstr( mlp, 108);
	if (_res) _res = _doevalstr( mlp, 109);
	if (_res) _res = _doevalstr( mlp, 110);
	if (_res) _res = MLPutSymbol( mlp, "End");
	if (_res) _res = MLFlush( mlp);
	return _res;
} /* MLInstall */


#if MLPROTOTYPES
int MLDoCallPacket( MLINK mlp)
#else
int MLDoCallPacket( mlp) MLINK mlp;
#endif
{
	return _MLDoCallPacket( mlp, _tramps, 111);
} /* MLDoCallPacket */

/******************************* begin trailer ********************************/

#ifndef EVALSTRS_AS_BYTESTRINGS
#	define EVALSTRS_AS_BYTESTRINGS 1
#endif

#if CARDOF_EVALSTRS
static int  _doevalstr( MLINK mlp, int n)
{
	long bytesleft, charsleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsnow;
#endif
	char **s, **p;
	char *t;

	s = (char **)evalstrs;
	while( n-- > 0){
		if( *s == 0) break;
		while( *s++ != 0){}
	}
	if( *s == 0) return 0;
	bytesleft = 0;
	charsleft = 0;
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft += bytesnow;
		charsleft += bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
		t = *p;
		charsleft -= MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
#endif
		++p;
	}


	MLPutNext( mlp, MLTKSTR);
#if EVALSTRS_AS_BYTESTRINGS
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		MLPut8BitCharacters( mlp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	MLPut7BitCount( mlp, charsleft, bytesleft);
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		t = *p;
		charsnow = bytesnow - MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
		charsleft -= charsnow;
		MLPut7BitCharacters(  mlp, charsleft, *p, bytesnow, charsnow);
		++p;
	}
#endif
	return MLError( mlp) == MLEOK;
}
#endif /* CARDOF_EVALSTRS */


static int  _definepattern( MLINK mlp, char* patt, char* args, int func_n)
{
	MLPutFunction( mlp, "DefineExternal", (long)3);
	  MLPutString( mlp, patt);
	  MLPutString( mlp, args);
	  MLPutInteger( mlp, func_n);
	return !MLError(mlp);
} /* _definepattern */


int _MLDoCallPacket( MLINK mlp, struct func functable[], int nfuncs)
{
#if MLINTERFACE >= 4
	int len;
#else
	long len;
#endif
	int n, res = 0;
	struct func* funcp;

	if( ! MLGetInteger( mlp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
#if MLINTERFACE >= 4
	&& ( ! MLTestHead(mlp, "List", &len)
#else
	&& ( ! MLCheckFunction(mlp, "List", &len)
#endif
	     || ( !funcp->manual && (len != funcp->f_nargs))
	     || (  funcp->manual && (len <  funcp->f_nargs))
	   )
	) goto L0;

	stdlink = mlp;
	res = (*funcp->f_func)( mlp);

L0:	if( res == 0)
		res = MLClearError( mlp) && MLPutSymbol( mlp, "$Failed");
	return res && MLEndPacket( mlp) && MLNewPacket( mlp);
} /* _MLDoCallPacket */


mlapi_packet MLAnswer( MLINK mlp)
{
	mlapi_packet pkt = 0;
#if MLINTERFACE >= 4
	int waitResult;

	while( ! MLDone && ! MLError(mlp)
		&& (waitResult = MLWaitForLinkActivity(mlp),waitResult) &&
		waitResult == MLWAITSUCCESS && (pkt = MLNextPacket(mlp), pkt) &&
		pkt == CALLPKT)
	{
		MLAbort = 0;
		if(! MLDoCallPacket(mlp))
			pkt = 0;
	}
#else
	while( !MLDone && !MLError(mlp) && (pkt = MLNextPacket(mlp), pkt) && pkt == CALLPKT){
		MLAbort = 0;
		if( !MLDoCallPacket(mlp)) pkt = 0;
	}
#endif
	MLAbort = 0;
	return pkt;
} /* MLAnswer */



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

static int refuse_to_be_a_frontend( MLINK mlp)
{
	int pkt;

	MLPutFunction( mlp, "EvaluatePacket", 1);
	  MLPutFunction( mlp, "Module", 2);
	    MLPutFunction( mlp, "List", 1);
		  MLPutFunction( mlp, "Set", 2);
		    MLPutSymbol( mlp, "me");
	        MLPutSymbol( mlp, "$ParentLink");
	  MLPutFunction( mlp, "CompoundExpression", 3);
	    MLPutFunction( mlp, "Set", 2);
	      MLPutSymbol( mlp, "$ParentLink");
	      MLTransferExpression( mlp, mlp);
	    MLPutFunction( mlp, "Message", 2);
	      MLPutFunction( mlp, "MessageName", 2);
	        MLPutSymbol( mlp, "$ParentLink");
	        MLPutString( mlp, "notfe");
	      MLPutSymbol( mlp, "me");
	    MLPutSymbol( mlp, "me");
	MLEndPacket( mlp);

	while( (pkt = MLNextPacket( mlp), pkt) && pkt != SUSPENDPKT)
		MLNewPacket( mlp);
	MLNewPacket( mlp);
	return MLError( mlp) == MLEOK;
}


int MLEvaluate( MLINK mlp, char* s)
{
	if( MLAbort) return 0;
	return MLPutFunction( mlp, "EvaluatePacket", 1L)
		&& MLPutFunction( mlp, "ToExpression", 1L)
		&& MLPutString( mlp, s)
		&& MLEndPacket( mlp);
} /* MLEvaluate */


int MLEvaluateString( MLINK mlp, char* s)
{
	int pkt;
	if( MLAbort) return 0;
	if( MLEvaluate( mlp, s)){
		while( (pkt = MLAnswer( mlp), pkt) && pkt != RETURNPKT)
			MLNewPacket( mlp);
		MLNewPacket( mlp);
	}
	return MLError( mlp) == MLEOK;
} /* MLEvaluateString */


#if MLINTERFACE >= 3
MLMDEFN( void, MLDefaultHandler, ( MLINK mlp, int message, int n))
#else
MLMDEFN( void, MLDefaultHandler, ( MLINK mlp, unsigned long message, unsigned long n))
#endif /* MLINTERFACE >= 3 */
{
	mlp = (MLINK)0; /* suppress unused warning */
	n = 0; /* suppress unused warning */

	switch (message){
	case MLTerminateMessage:
		MLDone = 1;
	case MLInterruptMessage:
	case MLAbortMessage:
		MLAbort = 1;
	default:
		return;
	}
}

#if MLINTERFACE >= 3
static int _MLMain( char **argv, char **argv_end, char *commandline)
#else
static int _MLMain( charpp_ct argv, charpp_ct argv_end, charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
{
	MLINK mlp;
#if MLINTERFACE >= 3
	int err;
#else
	long err;
#endif /* MLINTERFACE >= 3 */

#if (DARWIN_MATHLINK && CARBON_MPREP)
	if( !init_macintosh()) goto R0;
#endif /* (DARWIN_MATHLINK && CARBON_MPREP) */

#if MLINTERFACE >= 4
	if( !stdenv)
		stdenv = MLInitialize( (MLEnvironmentParameter)0);
#else
	if( !stdenv)
		stdenv = MLInitialize( (MLParametersPointer)0);
#endif

	if( stdenv == (MLEnvironment)0) goto R0;

#if (DARWIN_MATHLINK && CARBON_MPREP)
#if MLINTERFACE >= 3
	if( !stdyielder)
		stdyielder = (MLYieldFunctionObject)MLDefaultYielder;
#else
	if( !stdyielder)
		stdyielder = MLCreateYieldFunction( stdenv, NewMLYielderProc( MLDefaultYielder), 0);
#endif /* MLINTERFACE >= 3 */
#endif /* (DARWIN_MATHLINK && CARBON_MPREP)*/

#if MLINTERFACE >= 3
	if( !stdhandler)
		stdhandler = (MLMessageHandlerObject)MLDefaultHandler;
#else
	if( !stdhandler)
		stdhandler = MLCreateMessageHandler( stdenv, NewMLHandlerProc( MLDefaultHandler), 0);
#endif /* MLINTERFACE >= 3 */

#if (DARWIN_MATHLINK && CARBON_MPREP)
        MLSetDialogFunction(stdenv, MLRequestToInteractFunction, NewMLRequestToInteractProc(MLDontPermit_darwin));

	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
#if MLINTERFACE >= 3
			: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
#else
			: MLOpenArgv( stdenv, argv, argv_end, &err);
#endif

        MLSetDialogFunction(stdenv, MLRequestToInteractFunction, NewMLRequestToInteractProc(MLPermit_darwin));

	if( mlp == (MLINK)0){
                        mlp = commandline
                                ? MLOpenString( stdenv, commandline, &err)
#if MLINTERFACE < 3
                                : MLOpenArgv( stdenv, argv, argv_end, &err);
#else
                                : MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
#endif
        }
#else /* !(DARWIN_MATHLINK && CARBON_MPREP)*/
	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
#if MLINTERFACE < 3
		: MLOpenArgv( stdenv, argv, argv_end, &err);
#else
		: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
#endif
#endif /* (DARWIN_MATHLINK && CARBON_MPREP)*/

	if( mlp == (MLINK)0){
		MLAlert( stdenv, MLErrorString( stdenv, err));
		goto R1;
	}

	if( stdyielder) MLSetYieldFunction( mlp, stdyielder);
	if( stdhandler) MLSetMessageHandler( mlp, stdhandler);

	if( MLInstall( mlp))
		while( MLAnswer( mlp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( mlp)) break;
		}

	MLClose( mlp);
R1:	MLDeinitialize( stdenv);
	stdenv = (MLEnvironment)0;
R0:	return !MLDone;
} /* _MLMain */


#if MLINTERFACE >= 3
int MLMainString( char *commandline)
#else
int MLMainString( charp_ct commandline)
#endif /* MLINTERFACE >= 3 */
{
#if MLINTERFACE >= 3
	return _MLMain( (char **)0, (char **)0, commandline);
#else
	return _MLMain( (charpp_ct)0, (charpp_ct)0, commandline);
#endif /* MLINTERFACE >= 3 */
}

int MLMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
{   
	static char FAR * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
	return _MLMain( far_argv, far_argv + count, (charp_ct)0);

}

#if MLINTERFACE >= 3
int MLMain( int argc, char ** argv)
#else
int MLMain( int argc, charpp_ct argv)
#endif /* MLINTERFACE >= 3 */
{
#if MLINTERFACE >= 3
 	return _MLMain( argv, argv + argc, (char *)0);
#else
 	return _MLMain( argv, argv + argc, (charp_ct)0);
#endif /* MLINTERFACE >= 3 */
}
