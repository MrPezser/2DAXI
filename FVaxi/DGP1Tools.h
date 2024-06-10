//
// Created by Tsail on 5/29/2024.
//

#ifndef FVAXI_DGP1TOOLS_H
#define FVAXI_DGP1TOOLS_H

#include "Indexing.h"
#include "StateVariables.h"

void get_u_val(const double* unk, State var, Thermo air,const double* ux, const double* uy, double xsi, double eta, double* uout);
void get_u_val_standardrecon(const double* unk,const double* ux, const double* uy, double xsi, double eta, double* uout);
void DGP1_volume_integral(int nx, int ny, double vol, double* xfa, double* yfa, double* geoel, double* unk, State* ElemVar,
                          double* duxdt, double* duydt);
void DGP1_xsi_face_integral(int ieL, int ieR, int iuL, int iuR,double* unk, State* ElemVar, double* ux, double* uy,
                            const double* yCenter, Thermo air, double rFace, double* fNormal, double len,
                            double* rhsel, double* rhselx, double* rhsely);
void DGP1_eta_face_integral(int ieL, int ieR, int iuL, int iuR,double* unk, State* ElemVar, double* ux, double* uy,
                            const double* yCenter, Thermo air, double rFace, double* fNormal, double len,
                            double* rhsel, double* rhselx, double* rhsely);
void DGP1_boundary_face_integral(int ieIn, int ieEx, int iuIn, int iuEx,double* unk, State* ElemVar, double* ux, double* uy,
                                 int iFaceType, double* unkExt, State* EVExt, double yCenter, Thermo air, double rFace,
                                 double* fNormal, double* fNormalL, double* fNormalR, double len,
                                 double* rhsel, double* rhselx, double* rhsely);
void DGP1_ghost_cell_generator(int nx, int ny, double* unk, double* ux, double* uy, State* ElemVar, Thermo air, int* ibound,
                          double* geofa, double* uFS, double* uGBot, double* uGTop, double* uGLeft, double* uGRight,
                          State* BotVar, State* TopVar, State* LeftVar, State* RightVar);

#endif //FVAXI_DGP1TOOLS_H
