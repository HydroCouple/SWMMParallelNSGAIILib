/*! \file   swmmparallelnsgaii.h
 *  \author Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version   1.0.0.0
 *  \section   Description
 *  This header files contains is part of the swmmparallelnsgaii library.
 *  It defines the macro used for exporting externally callable functions in the shared library.
 *  \section License
 *  swmmparallelnsgaii.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  The swmmparallelnsgaii library and its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \warning
 *  \todo
 */


#ifndef SWMMPARALLELNSGAII_GLOBAL_H
#define SWMMPARALLELNSGAII_GLOBAL_H


#  ifdef Q_OS_WIN
#    define SWMMParallelNSGAII_EXPORT     __declspec(dllexport)
#    define SWMMParallelNSGAII_IMPORT     __declspec(dllimport)
#  else
#    define SWMMParallelNSGAII_EXPORT     __attribute__((visibility("default")))
#    define SWMMParallelNSGAII_IMPORT     __attribute__((visibility("default")))
#    define SWMMParallelNSGAII_HIDDEN     __attribute__((visibility("hidden")))
#  endif

#endif // SWMMPARALLELNSGAII_GLOBAL_H
