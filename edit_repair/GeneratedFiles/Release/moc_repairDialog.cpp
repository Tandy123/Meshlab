/****************************************************************************
** Meta object code from reading C++ file 'repairDialog.h'
**
** Created: Tue Mar 7 14:48:48 2017
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../repairDialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'repairDialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_RepairDialog[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      14,   13,   13,   13, 0x05,
      24,   13,   13,   13, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_RepairDialog[] = {
    "RepairDialog\0\0closing()\0"
    "updateMeshSetVisibilities()\0"
};

const QMetaObject RepairDialog::staticMetaObject = {
    { &QDockWidget::staticMetaObject, qt_meta_stringdata_RepairDialog,
      qt_meta_data_RepairDialog, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &RepairDialog::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *RepairDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *RepairDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_RepairDialog))
        return static_cast<void*>(const_cast< RepairDialog*>(this));
    return QDockWidget::qt_metacast(_clname);
}

int RepairDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDockWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: closing(); break;
        case 1: updateMeshSetVisibilities(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void RepairDialog::closing()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void RepairDialog::updateMeshSetVisibilities()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
