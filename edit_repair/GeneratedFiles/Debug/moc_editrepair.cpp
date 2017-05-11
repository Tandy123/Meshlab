/****************************************************************************
** Meta object code from reading C++ file 'editrepair.h'
**
** Created: Sun Feb 5 21:51:09 2017
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../editrepair.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'editrepair.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_EditRepairPlugin[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      39,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      18,   17,   17,   17, 0x0a,
      32,   17,   17,   17, 0x0a,
      53,   17,   17,   17, 0x0a,
      75,   17,   17,   17, 0x0a,
      97,   17,   17,   17, 0x0a,
     121,   17,   17,   17, 0x0a,
     143,   17,   17,   17, 0x0a,
     168,   17,   17,   17, 0x0a,
     193,   17,   17,   17, 0x0a,
     218,   17,   17,   17, 0x0a,
     241,   17,   17,   17, 0x0a,
     257,   17,   17,   17, 0x0a,
     273,   17,   17,   17, 0x0a,
     287,   17,   17,   17, 0x0a,
     302,   17,   17,   17, 0x0a,
     321,   17,   17,   17, 0x0a,
     338,   17,   17,   17, 0x0a,
     356,   17,   17,   17, 0x0a,
     371,   17,   17,   17, 0x0a,
     390,   17,   17,   17, 0x0a,
     405,   17,   17,   17, 0x0a,
     425,   17,   17,   17, 0x0a,
     445,   17,   17,   17, 0x0a,
     465,   17,   17,   17, 0x0a,
     489,   17,   17,   17, 0x0a,
     506,   17,   17,   17, 0x0a,
     530,   17,   17,   17, 0x0a,
     555,   17,   17,   17, 0x0a,
     580,   17,   17,   17, 0x0a,
     605,   17,   17,   17, 0x0a,
     622,   17,   17,   17, 0x0a,
     649,   17,   17,   17, 0x0a,
     670,   17,   17,   17, 0x0a,
     687,   17,   17,   17, 0x0a,
     703,   17,   17,   17, 0x0a,
     720,   17,   17,   17, 0x0a,
     740,   17,   17,   17, 0x0a,
     758,   17,   17,   17, 0x0a,
     768,   17,   17,   17, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_EditRepairPlugin[] = {
    "EditRepairPlugin\0\0slotProject()\0"
    "slotDeformationPro()\0slotDeformationPro2()\0"
    "slotDeformationPro3()\0slotDeformationNonPro()\0"
    "slotDeformationAuto()\0slotDeformationNonPro2()\0"
    "slotDeformationNonPro3()\0"
    "slotDeformationNonPro4()\0"
    "slotDeformationAuto2()\0slotDeletePro()\0"
    "slotImportSet()\0slotFillGap()\0"
    "slotFillGap2()\0slotFindBoundary()\0"
    "slotSortBorder()\0slotSortBorder2()\0"
    "slotCopyMesh()\0slotImportSetTem()\0"
    "slotAlignTem()\0slotAdjustX(double)\0"
    "slotAdjustY(double)\0slotAdjustZ(double)\0"
    "slotAdjustScale(double)\0slotSaveAdjust()\0"
    "slotProjectBarycentre()\0"
    "slotProjectBarycentre2()\0"
    "slotProjectBarycentre3()\0"
    "slotProjectBarycentre4()\0slotExportMesh()\0"
    "slotDeformationLaplacian()\0"
    "slotGetBoundaryTem()\0slotSmoothRoot()\0"
    "slotSmoothGap()\0slotSmoothGap2()\0"
    "slotDefomationTem()\0slotRepairByTem()\0"
    "slotCut()\0slotMerge()\0"
};

const QMetaObject EditRepairPlugin::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_EditRepairPlugin,
      qt_meta_data_EditRepairPlugin, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &EditRepairPlugin::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *EditRepairPlugin::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *EditRepairPlugin::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_EditRepairPlugin))
        return static_cast<void*>(const_cast< EditRepairPlugin*>(this));
    if (!strcmp(_clname, "MeshEditInterface"))
        return static_cast< MeshEditInterface*>(const_cast< EditRepairPlugin*>(this));
    if (!strcmp(_clname, "vcg.meshlab.MeshEditInterface/1.0"))
        return static_cast< MeshEditInterface*>(const_cast< EditRepairPlugin*>(this));
    return QObject::qt_metacast(_clname);
}

int EditRepairPlugin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: slotProject(); break;
        case 1: slotDeformationPro(); break;
        case 2: slotDeformationPro2(); break;
        case 3: slotDeformationPro3(); break;
        case 4: slotDeformationNonPro(); break;
        case 5: slotDeformationAuto(); break;
        case 6: slotDeformationNonPro2(); break;
        case 7: slotDeformationNonPro3(); break;
        case 8: slotDeformationNonPro4(); break;
        case 9: slotDeformationAuto2(); break;
        case 10: slotDeletePro(); break;
        case 11: slotImportSet(); break;
        case 12: slotFillGap(); break;
        case 13: slotFillGap2(); break;
        case 14: slotFindBoundary(); break;
        case 15: slotSortBorder(); break;
        case 16: slotSortBorder2(); break;
        case 17: slotCopyMesh(); break;
        case 18: slotImportSetTem(); break;
        case 19: slotAlignTem(); break;
        case 20: slotAdjustX((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 21: slotAdjustY((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 22: slotAdjustZ((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 23: slotAdjustScale((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 24: slotSaveAdjust(); break;
        case 25: slotProjectBarycentre(); break;
        case 26: slotProjectBarycentre2(); break;
        case 27: slotProjectBarycentre3(); break;
        case 28: slotProjectBarycentre4(); break;
        case 29: slotExportMesh(); break;
        case 30: slotDeformationLaplacian(); break;
        case 31: slotGetBoundaryTem(); break;
        case 32: slotSmoothRoot(); break;
        case 33: slotSmoothGap(); break;
        case 34: slotSmoothGap2(); break;
        case 35: slotDefomationTem(); break;
        case 36: slotRepairByTem(); break;
        case 37: slotCut(); break;
        case 38: slotMerge(); break;
        default: ;
        }
        _id -= 39;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
