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
#include <SofaMJEDFEM/config.h>
#include <sofa/core/ObjectFactory.h>
using sofa::core::ObjectFactory;

namespace sofa::component
{

extern "C" {
    SOFA_MJED_FEM_API void initExternalModule();
    SOFA_MJED_FEM_API const char* getModuleName();
    SOFA_MJED_FEM_API const char* getModuleVersion();
    SOFA_MJED_FEM_API const char* getModuleLicense();
    SOFA_MJED_FEM_API const char* getModuleDescription();
    SOFA_MJED_FEM_API const char* getModuleComponentList();
}

void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleName()
{
    return sofa_tostring(SOFA_TARGET);
}

const char* getModuleVersion()
{
    return sofa_tostring(SOFAMJEDFEM_VERSION);
}

const char* getModuleLicense()
{
    return "LGPL";
}

const char* getModuleDescription()
{
    return "This plugin implements the Multiplicative Jacobian Energy Decomposition or MJED method.";
}

const char* getModuleComponentList()
{
    /// string containing the names of the classes provided by the plugin
    static std::string classes = ObjectFactory::getInstance()->listClassesFromTarget(sofa_tostring(SOFA_TARGET));
    return classes.c_str();
}

} // namespace sofa::component
