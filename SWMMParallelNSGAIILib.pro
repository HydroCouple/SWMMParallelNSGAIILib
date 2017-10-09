VERSION = 1.0
TARGET = SWMMParallelNSGAIILib
TEMPLATE = lib
QT       -= gui
CONFIG += c++11 -stdlib=libstdc++


#DEFINES += USE_OPENMP

INCLUDEPATH += .\
               ./include \
               ./include/swmm_core

HEADERS += \
    include/swmmparallelnsgaii_global.h \
    include/swmmparallelnsgaii.h \
    include/swmm_core/headers.h \
    include/swmm_core/consts.h \
    include/swmm_core/datetime.h \
    include/swmm_core/enums.h \
    include/swmm_core/error.h \
    include/swmm_core/exfil.h \
    include/swmm_core/findroot.h \
    include/swmm_core/funcs.h \
    include/swmm_core/globals.h \
    include/swmm_core/hash.h \
    include/swmm_core/infil.h \
    include/swmm_core/keywords.h \
    include/swmm_core/lid.h \
    include/swmm_core/macros.h \
    include/swmm_core/mathexpr.h \
    include/swmm_core/mempool.h \
    include/swmm_core/objects.h \
    include/swmm_core/odesolve.h \
    include/swmm_core/swmm5.h \
    include/swmm_core/text.h \
    include/swmm_core/swmm5_iface.h

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
    src/swmmparallelnsgaii.cpp \
    src/swmm_core/iface.c \
    src/swmm_core/error.c \
    src/swmm_core/odesolve.c \
    src/swmm_core/controls.c \
    src/swmm_core/inflow.c \
    src/swmm_core/exfil.c \
    src/swmm_core/snow.c \
    src/swmm_core/datetime.c \
    src/swmm_core/runoff.c \
    src/swmm_core/forcmain.c \
    src/swmm_core/gage.c \
    src/swmm_core/lid.c \
    src/swmm_core/massbal.c \
    src/swmm_core/infil.c \
    src/swmm_core/climate.c \
    src/swmm_core/hash.c \
    src/swmm_core/dwflow.c \
    src/swmm_core/culvert.c \
    src/swmm_core/keywords.c \
    src/swmm_core/hotstart.c \
    src/swmm_core/gwater.c \
    src/swmm_core/table.c \
    src/swmm_core/kinwave.c \
    src/swmm_core/lidproc.c \
    src/swmm_core/swmm5.c \
    src/swmm_core/findroot.c \
    src/swmm_core/node.c \
    src/swmm_core/rain.c \
    src/swmm_core/output.c \
    src/swmm_core/inputrpt.c \
    src/swmm_core/project.c \
    src/swmm_core/dynwave.c \
    src/swmm_core/qualrout.c \
    src/swmm_core/shape.c \
    src/swmm_core/rdii.c \
    src/swmm_core/mempool.c \
    src/swmm_core/routing.c \
    src/swmm_core/input.c \
    src/swmm_core/mathexpr.c \
    src/swmm_core/stats.c \
    src/swmm_core/report.c \
    src/swmm_core/landuse.c \
    src/swmm_core/surfqual.c \
    src/swmm_core/subcatch.c \
    src/swmm_core/toposort.c \
    src/swmm_core/link.c \
    src/swmm_core/transect.c \
    src/swmm_core/statsrpt.c \
    src/swmm_core/treatmnt.c \
    src/swmm_core/flowrout.c \
    src/swmm_core/roadway.c \
    src/swmm_core/xsect.c \
    src/swmm_core/swmm5_iface.c

DISTFILES += \
    include/swmm_core/xsect.dat


