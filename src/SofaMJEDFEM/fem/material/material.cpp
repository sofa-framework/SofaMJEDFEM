/******************************************************************************
*                             SOFA MJED FEM plugin                            *
*                                (c) 2006 Inria                               *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: see Authors.md                                                     *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//used only to test compilation and dependencies in SOFA. (TODO: see if it has to be removed, or to be removed as soon as another .cpp include them)

#include <SofaMiscFem/BoyceAndArruda.h>
#include <SofaMiscFem/MooneyRivlin.h>
#include <SofaMiscFem/STVenantKirchhoff.h>
#include <SofaMiscFem/HyperelasticMaterial.h>
#include <SofaMiscFem/VerondaWestman.h>
#include <SofaMiscFem/Costa.h>
#include <SofaMJEDFEM/fem/material/NeoHookean.h>
#include <SofaMJEDFEM/fem/material/NeoHookeanIsotropicMJED.h>
#include <SofaMJEDFEM/fem/material/Ogden.h>
#include <SofaMJEDFEM/fem/material/BoyceAndArrudaMJED.h>
#include <SofaMJEDFEM/fem/material/MooneyRivlinMJED.h>
#include <SofaMJEDFEM/fem/material/STVenantKirchhoffMJED.h>
#include <SofaMJEDFEM/fem/material/HyperelasticMaterialMJED.h>
#include <SofaMJEDFEM/fem/material/NeoHookeanMJED.h>
#include <SofaMJEDFEM/fem/material/VerondaWestmanMJED.h>
#include <SofaMJEDFEM/fem/material/CostaMJED.h>
#include <SofaMJEDFEM/fem/material/OgdenMJED.h>

