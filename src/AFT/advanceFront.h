/* 
 * File:   advanceFront.h
 * Author: orhan
 *
 * Created on July 18, 2015, 10:51 PM
 */

#ifndef ADVANCEFRONT_H
#define	ADVANCEFRONT_H

#define tolPreferExistingPoint 0.0
double scoreAB;
double scoreNP1;
double scoreNP2;
double score;
bool A_CPX_exists;
bool B_CPX_exists;
int iA_CPX;
int iB_CPX;
bool elAB;
bool elNP1;
bool elNP2;
Point cnp1;
Point cnp2;
int iCnp1;
int iCnp2;
int iChosenPoint;
bool isNewPoint;
//bool A_CPX_exists_AB;
//bool B_CPX_exists_AB;
//int iA_CPX_AB;
//int iB_CPX_AB;
bool A_CPX_exists_NP1;
bool B_CPX_exists_NP1;
int iA_CPX_NP1;
int iB_CPX_NP1;
bool A_CPX_exists_NP2;
bool B_CPX_exists_NP2;
int iA_CPX_NP2;
int iB_CPX_NP2;
double pdis;

#endif	/* ADVANCEFRONT_H */

