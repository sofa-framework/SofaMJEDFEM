/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_COMPONENT_FORCEFIELD_MJEDTETRAHEDRALFORCEFIELD_CPP

#include <SofaMJEDFEM/MJEDTetrahedralForceField.inl>

#include <SofaMJEDFEM/config.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>

#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>

#include <string.h>
#include <iostream>

namespace sofa::component::forcefield
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(MJEDTetrahedralForceField)

// Register in the Factory
int MJEDTetrahedralForceFieldClass = core::RegisterObject("Generic Tetrahedral finite elements")
#ifndef SOFA_FLOAT
.add< MJEDTetrahedralForceField<sofa::defaulttype::Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< MJEDTetrahedralForceField<Vec3fTypes> >()
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_MJED_FEM_API MJEDTetrahedralForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_MJED_FEM_API MJEDTetrahedralForceField<Vec3fTypes>;
#endif

} // namespace sofa::component::forcefield
