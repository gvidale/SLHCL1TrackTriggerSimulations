/*

 * functions.h
 *
 *  Created on: 08 giu 2017
 *      Author: gvidale
 */

#ifndef MYMACRO_FUNCTIONS_H_
#define MYMACRO_FUNCTIONS_H_

#include <TROOT.h>

Int_t* phi_zeta(Int_t modId){

	Int_t* lpz= new Int_t[3];
//	Int_t layer, phi, zeta;

	lpz[0] = modId/10000;
	lpz[1] = ((modId)/100)%100;
	lpz[2] = (modId)%100;

	return lpz;
}





#endif /* MYMACRO_FUNCTIONS_H_ */
