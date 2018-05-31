#Author Caleb Amoa Buahin
#Email caleb.buahin@gmail.com
#Date 2016
#License GNU General Public License (see <http://www.gnu.org/licenses/> for details).

VERSION = 1.0
TARGET = SWMMParallelNSGAIILib
TEMPLATE = lib
QT       -= gui
CONFIG += c++11 -stdlib=libstdc++


#DEFINES += USE_OPENMP

INCLUDEPATH += .\
               ./include \
               ../SWMM/5.1.012/include\
               ../SWMM/swmm5_iface/include

HEADERS += \
    include/swmmparallelnsgaii_global.h \
    include/swmmparallelnsgaii.h

CONFIG(debug, debug|release){
   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui
}

CONFIG(release, debug|release){
    DESTDIR = lib
    RELEASE_EXTRAS = ./build/release
    OBJECTS_DIR = $$RELEASE_EXTRAS/.obj
    MOC_DIR = $$RELEASE_EXTRAS/.moc
    RCC_DIR = $$RELEASE_EXTRAS/.qrc
    UI_DIR = $$RELEASE_EXTRAS/.ui
}

SOURCES += \
    src/swmmparallelnsgaii.cpp


