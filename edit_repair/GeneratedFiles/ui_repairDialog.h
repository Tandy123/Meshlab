/********************************************************************************
** Form generated from reading UI file 'repairDialog.ui'
**
** Created: Tue Mar 7 14:48:50 2017
**      by: Qt User Interface Compiler version 4.7.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_REPAIRDIALOG_H
#define UI_REPAIRDIALOG_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFormLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Form
{
public:
    QPushButton *projectButton;
    QPushButton *deformatinoNonProButton;
    QPushButton *alignTemButton;
    QPushButton *importSetTemButton;
    QDoubleSpinBox *zAdjustdoubleSpinBox;
    QPushButton *saveAdjustButton;
    QDoubleSpinBox *xAdjustdoubleSpinBox;
    QDoubleSpinBox *yAdjustdoubleSpinBox;
    QDoubleSpinBox *scaledoubleSpinBox;
    QPushButton *exportMeshButton;
    QPushButton *importSetButton_2;
    QPushButton *defomationTemButton;
    QPushButton *repairButton;
    QPushButton *deformationButton;
    QPushButton *deformationAutoButton;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout;
    QPushButton *importSetButton;
    QPushButton *projectBarycentreButton;
    QPushButton *deformationButton2;
    QHBoxLayout *horizontalLayout;
    QPushButton *getBoundaryTemButton;
    QCheckBox *checkBoxShowPoints;
    QPushButton *deformationAutoButton2;
    QPushButton *smoothRootButton;
    QPushButton *deleteProButton;
    QPushButton *copyMeshButton;
    QPushButton *findBoundaryButton;
    QPushButton *sortBorderButton;
    QPushButton *fillGapButton;
    QPushButton *smoothGapButton;
    QPushButton *cutButton;
    QPushButton *mergeButton;
    QPushButton *findNeastestPointsButton;
    QPushButton *normalizeButton;
    QLabel *label;
    QCheckBox *checkBoxHideAxis;
    QDoubleSpinBox *doubleSpinBoxPerStep;
    QPushButton *exportPButton;
    QSpinBox *spinBoxCutTimes;
    QLabel *label_2;
    QDoubleSpinBox *doubleSpinBoxPlaneOffset;
    QLabel *label_3;
    QLabel *toothinfo;
    QPushButton *exportPButton2;
    QPushButton *errorMapButton;
    QLabel *label_8;
    QWidget *layoutWidget1;
    QFormLayout *formLayout;
    QLabel *label_4;
    QDoubleSpinBox *doubleSpinBox;
    QLabel *label_5;
    QDoubleSpinBox *doubleSpinBox_2;
    QLabel *label_6;
    QDoubleSpinBox *doubleSpinBox_3;
    QLabel *label_7;
    QDoubleSpinBox *doubleSpinBox_4;
    QLabel *errorinfo;
    QLabel *errorinfo_2;
    QPushButton *invertNormalButton;
    QPushButton *invertNormalFaceButton;
    QPushButton *extendBorderButton;
    QPushButton *linkPointsButton;
    QPushButton *importInsectionPButton;

    void setupUi(QWidget *Form)
    {
        if (Form->objectName().isEmpty())
            Form->setObjectName(QString::fromUtf8("Form"));
        Form->resize(317, 601);
        projectButton = new QPushButton(Form);
        projectButton->setObjectName(QString::fromUtf8("projectButton"));
        projectButton->setGeometry(QRect(210, 730, 151, 23));
        deformatinoNonProButton = new QPushButton(Form);
        deformatinoNonProButton->setObjectName(QString::fromUtf8("deformatinoNonProButton"));
        deformatinoNonProButton->setGeometry(QRect(200, 760, 161, 23));
        alignTemButton = new QPushButton(Form);
        alignTemButton->setObjectName(QString::fromUtf8("alignTemButton"));
        alignTemButton->setGeometry(QRect(410, 640, 111, 23));
        importSetTemButton = new QPushButton(Form);
        importSetTemButton->setObjectName(QString::fromUtf8("importSetTemButton"));
        importSetTemButton->setGeometry(QRect(410, 610, 111, 23));
        zAdjustdoubleSpinBox = new QDoubleSpinBox(Form);
        zAdjustdoubleSpinBox->setObjectName(QString::fromUtf8("zAdjustdoubleSpinBox"));
        zAdjustdoubleSpinBox->setGeometry(QRect(39, 735, 62, 22));
        zAdjustdoubleSpinBox->setMinimum(-99.99);
        zAdjustdoubleSpinBox->setSingleStep(0.2);
        saveAdjustButton = new QPushButton(Form);
        saveAdjustButton->setObjectName(QString::fromUtf8("saveAdjustButton"));
        saveAdjustButton->setGeometry(QRect(30, 770, 75, 23));
        xAdjustdoubleSpinBox = new QDoubleSpinBox(Form);
        xAdjustdoubleSpinBox->setObjectName(QString::fromUtf8("xAdjustdoubleSpinBox"));
        xAdjustdoubleSpinBox->setGeometry(QRect(39, 675, 62, 22));
        xAdjustdoubleSpinBox->setMinimum(-99.99);
        xAdjustdoubleSpinBox->setSingleStep(0.2);
        yAdjustdoubleSpinBox = new QDoubleSpinBox(Form);
        yAdjustdoubleSpinBox->setObjectName(QString::fromUtf8("yAdjustdoubleSpinBox"));
        yAdjustdoubleSpinBox->setGeometry(QRect(39, 705, 62, 22));
        yAdjustdoubleSpinBox->setMinimum(-99.99);
        yAdjustdoubleSpinBox->setSingleStep(0.2);
        scaledoubleSpinBox = new QDoubleSpinBox(Form);
        scaledoubleSpinBox->setObjectName(QString::fromUtf8("scaledoubleSpinBox"));
        scaledoubleSpinBox->setGeometry(QRect(110, 680, 62, 22));
        scaledoubleSpinBox->setMaximum(2);
        scaledoubleSpinBox->setSingleStep(0.05);
        scaledoubleSpinBox->setValue(1);
        exportMeshButton = new QPushButton(Form);
        exportMeshButton->setObjectName(QString::fromUtf8("exportMeshButton"));
        exportMeshButton->setGeometry(QRect(30, 800, 75, 23));
        importSetButton_2 = new QPushButton(Form);
        importSetButton_2->setObjectName(QString::fromUtf8("importSetButton_2"));
        importSetButton_2->setGeometry(QRect(150, 10, 151, 23));
        defomationTemButton = new QPushButton(Form);
        defomationTemButton->setObjectName(QString::fromUtf8("defomationTemButton"));
        defomationTemButton->setGeometry(QRect(410, 550, 101, 23));
        repairButton = new QPushButton(Form);
        repairButton->setObjectName(QString::fromUtf8("repairButton"));
        repairButton->setGeometry(QRect(421, 575, 81, 23));
        deformationButton = new QPushButton(Form);
        deformationButton->setObjectName(QString::fromUtf8("deformationButton"));
        deformationButton->setGeometry(QRect(210, 670, 151, 23));
        deformationAutoButton = new QPushButton(Form);
        deformationAutoButton->setObjectName(QString::fromUtf8("deformationAutoButton"));
        deformationAutoButton->setGeometry(QRect(210, 700, 151, 23));
        layoutWidget = new QWidget(Form);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(400, 10, 239, 511));
        verticalLayout = new QVBoxLayout(layoutWidget);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        importSetButton = new QPushButton(layoutWidget);
        importSetButton->setObjectName(QString::fromUtf8("importSetButton"));

        verticalLayout->addWidget(importSetButton);

        projectBarycentreButton = new QPushButton(layoutWidget);
        projectBarycentreButton->setObjectName(QString::fromUtf8("projectBarycentreButton"));

        verticalLayout->addWidget(projectBarycentreButton);

        deformationButton2 = new QPushButton(layoutWidget);
        deformationButton2->setObjectName(QString::fromUtf8("deformationButton2"));

        verticalLayout->addWidget(deformationButton2);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        getBoundaryTemButton = new QPushButton(layoutWidget);
        getBoundaryTemButton->setObjectName(QString::fromUtf8("getBoundaryTemButton"));

        horizontalLayout->addWidget(getBoundaryTemButton);

        checkBoxShowPoints = new QCheckBox(layoutWidget);
        checkBoxShowPoints->setObjectName(QString::fromUtf8("checkBoxShowPoints"));

        horizontalLayout->addWidget(checkBoxShowPoints);


        verticalLayout->addLayout(horizontalLayout);

        deformationAutoButton2 = new QPushButton(layoutWidget);
        deformationAutoButton2->setObjectName(QString::fromUtf8("deformationAutoButton2"));

        verticalLayout->addWidget(deformationAutoButton2);

        smoothRootButton = new QPushButton(layoutWidget);
        smoothRootButton->setObjectName(QString::fromUtf8("smoothRootButton"));

        verticalLayout->addWidget(smoothRootButton);

        deleteProButton = new QPushButton(layoutWidget);
        deleteProButton->setObjectName(QString::fromUtf8("deleteProButton"));

        verticalLayout->addWidget(deleteProButton);

        copyMeshButton = new QPushButton(layoutWidget);
        copyMeshButton->setObjectName(QString::fromUtf8("copyMeshButton"));

        verticalLayout->addWidget(copyMeshButton);

        findBoundaryButton = new QPushButton(layoutWidget);
        findBoundaryButton->setObjectName(QString::fromUtf8("findBoundaryButton"));

        verticalLayout->addWidget(findBoundaryButton);

        sortBorderButton = new QPushButton(layoutWidget);
        sortBorderButton->setObjectName(QString::fromUtf8("sortBorderButton"));

        verticalLayout->addWidget(sortBorderButton);

        fillGapButton = new QPushButton(layoutWidget);
        fillGapButton->setObjectName(QString::fromUtf8("fillGapButton"));

        verticalLayout->addWidget(fillGapButton);

        smoothGapButton = new QPushButton(layoutWidget);
        smoothGapButton->setObjectName(QString::fromUtf8("smoothGapButton"));

        verticalLayout->addWidget(smoothGapButton);

        cutButton = new QPushButton(Form);
        cutButton->setObjectName(QString::fromUtf8("cutButton"));
        cutButton->setGeometry(QRect(10, 110, 291, 23));
        mergeButton = new QPushButton(Form);
        mergeButton->setObjectName(QString::fromUtf8("mergeButton"));
        mergeButton->setGeometry(QRect(0, 610, 112, 34));
        findNeastestPointsButton = new QPushButton(Form);
        findNeastestPointsButton->setObjectName(QString::fromUtf8("findNeastestPointsButton"));
        findNeastestPointsButton->setGeometry(QRect(10, 520, 141, 34));
        normalizeButton = new QPushButton(Form);
        normalizeButton->setObjectName(QString::fromUtf8("normalizeButton"));
        normalizeButton->setGeometry(QRect(120, 610, 112, 34));
        label = new QLabel(Form);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 10, 54, 21));
        checkBoxHideAxis = new QCheckBox(Form);
        checkBoxHideAxis->setObjectName(QString::fromUtf8("checkBoxHideAxis"));
        checkBoxHideAxis->setGeometry(QRect(150, 40, 121, 21));
        doubleSpinBoxPerStep = new QDoubleSpinBox(Form);
        doubleSpinBoxPerStep->setObjectName(QString::fromUtf8("doubleSpinBoxPerStep"));
        doubleSpinBoxPerStep->setGeometry(QRect(71, 10, 71, 22));
        doubleSpinBoxPerStep->setMinimum(-99);
        doubleSpinBoxPerStep->setSingleStep(0.1);
        doubleSpinBoxPerStep->setValue(-0.2);
        exportPButton = new QPushButton(Form);
        exportPButton->setObjectName(QString::fromUtf8("exportPButton"));
        exportPButton->setGeometry(QRect(150, 70, 151, 23));
        spinBoxCutTimes = new QSpinBox(Form);
        spinBoxCutTimes->setObjectName(QString::fromUtf8("spinBoxCutTimes"));
        spinBoxCutTimes->setGeometry(QRect(70, 70, 71, 22));
        spinBoxCutTimes->setMaximum(999);
        spinBoxCutTimes->setValue(60);
        label_2 = new QLabel(Form);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(10, 40, 54, 21));
        doubleSpinBoxPlaneOffset = new QDoubleSpinBox(Form);
        doubleSpinBoxPlaneOffset->setObjectName(QString::fromUtf8("doubleSpinBoxPlaneOffset"));
        doubleSpinBoxPlaneOffset->setGeometry(QRect(71, 40, 71, 22));
        doubleSpinBoxPlaneOffset->setMinimum(-99);
        doubleSpinBoxPlaneOffset->setValue(6);
        label_3 = new QLabel(Form);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(10, 70, 54, 21));
        toothinfo = new QLabel(Form);
        toothinfo->setObjectName(QString::fromUtf8("toothinfo"));
        toothinfo->setGeometry(QRect(10, 140, 291, 21));
        exportPButton2 = new QPushButton(Form);
        exportPButton2->setObjectName(QString::fromUtf8("exportPButton2"));
        exportPButton2->setGeometry(QRect(0, 640, 112, 34));
        errorMapButton = new QPushButton(Form);
        errorMapButton->setObjectName(QString::fromUtf8("errorMapButton"));
        errorMapButton->setGeometry(QRect(10, 330, 291, 34));
        label_8 = new QLabel(Form);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(10, 310, 161, 18));
        layoutWidget1 = new QWidget(Form);
        layoutWidget1->setObjectName(QString::fromUtf8("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(10, 170, 291, 129));
        formLayout = new QFormLayout(layoutWidget1);
        formLayout->setObjectName(QString::fromUtf8("formLayout"));
        formLayout->setContentsMargins(0, 0, 0, 0);
        label_4 = new QLabel(layoutWidget1);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        formLayout->setWidget(0, QFormLayout::LabelRole, label_4);

        doubleSpinBox = new QDoubleSpinBox(layoutWidget1);
        doubleSpinBox->setObjectName(QString::fromUtf8("doubleSpinBox"));
        doubleSpinBox->setEnabled(false);
        doubleSpinBox->setSingleStep(0.01);
        doubleSpinBox->setValue(0.04);

        formLayout->setWidget(0, QFormLayout::FieldRole, doubleSpinBox);

        label_5 = new QLabel(layoutWidget1);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        formLayout->setWidget(1, QFormLayout::LabelRole, label_5);

        doubleSpinBox_2 = new QDoubleSpinBox(layoutWidget1);
        doubleSpinBox_2->setObjectName(QString::fromUtf8("doubleSpinBox_2"));
        doubleSpinBox_2->setEnabled(false);
        doubleSpinBox_2->setSingleStep(0.01);
        doubleSpinBox_2->setValue(0.08);

        formLayout->setWidget(1, QFormLayout::FieldRole, doubleSpinBox_2);

        label_6 = new QLabel(layoutWidget1);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        formLayout->setWidget(2, QFormLayout::LabelRole, label_6);

        doubleSpinBox_3 = new QDoubleSpinBox(layoutWidget1);
        doubleSpinBox_3->setObjectName(QString::fromUtf8("doubleSpinBox_3"));
        doubleSpinBox_3->setEnabled(false);
        doubleSpinBox_3->setSingleStep(0.01);
        doubleSpinBox_3->setValue(0.12);

        formLayout->setWidget(2, QFormLayout::FieldRole, doubleSpinBox_3);

        label_7 = new QLabel(layoutWidget1);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        formLayout->setWidget(3, QFormLayout::LabelRole, label_7);

        doubleSpinBox_4 = new QDoubleSpinBox(layoutWidget1);
        doubleSpinBox_4->setObjectName(QString::fromUtf8("doubleSpinBox_4"));
        doubleSpinBox_4->setSingleStep(0.01);
        doubleSpinBox_4->setValue(0.16);

        formLayout->setWidget(3, QFormLayout::FieldRole, doubleSpinBox_4);

        errorinfo = new QLabel(Form);
        errorinfo->setObjectName(QString::fromUtf8("errorinfo"));
        errorinfo->setGeometry(QRect(10, 380, 291, 41));
        errorinfo_2 = new QLabel(Form);
        errorinfo_2->setObjectName(QString::fromUtf8("errorinfo_2"));
        errorinfo_2->setGeometry(QRect(10, 420, 291, 51));
        invertNormalButton = new QPushButton(Form);
        invertNormalButton->setObjectName(QString::fromUtf8("invertNormalButton"));
        invertNormalButton->setGeometry(QRect(250, 640, 121, 34));
        invertNormalFaceButton = new QPushButton(Form);
        invertNormalFaceButton->setObjectName(QString::fromUtf8("invertNormalFaceButton"));
        invertNormalFaceButton->setGeometry(QRect(120, 640, 121, 34));
        extendBorderButton = new QPushButton(Form);
        extendBorderButton->setObjectName(QString::fromUtf8("extendBorderButton"));
        extendBorderButton->setGeometry(QRect(161, 520, 141, 34));
        linkPointsButton = new QPushButton(Form);
        linkPointsButton->setObjectName(QString::fromUtf8("linkPointsButton"));
        linkPointsButton->setGeometry(QRect(10, 560, 291, 34));
        importInsectionPButton = new QPushButton(Form);
        importInsectionPButton->setObjectName(QString::fromUtf8("importInsectionPButton"));
        importInsectionPButton->setGeometry(QRect(10, 480, 291, 34));

        retranslateUi(Form);

        QMetaObject::connectSlotsByName(Form);
    } // setupUi

    void retranslateUi(QWidget *Form)
    {
        Form->setWindowTitle(QApplication::translate("Form", "Form", 0, QApplication::UnicodeUTF8));
        projectButton->setText(QApplication::translate("Form", "\346\212\225\345\275\261(\346\262\277\346\263\225\345\220\221)", 0, QApplication::UnicodeUTF8));
        deformatinoNonProButton->setText(QApplication::translate("Form", "\345\275\242\345\217\230\357\274\210\351\235\236\346\212\225\345\275\261\357\274\211", 0, QApplication::UnicodeUTF8));
        alignTemButton->setText(QApplication::translate("Form", "\347\247\273\345\212\250\346\250\241\345\236\213", 0, QApplication::UnicodeUTF8));
        importSetTemButton->setText(QApplication::translate("Form", "\345\257\274\345\205\245\346\250\241\346\235\277\345\261\200\351\203\250\345\235\220\346\240\207", 0, QApplication::UnicodeUTF8));
        saveAdjustButton->setText(QApplication::translate("Form", "\344\277\235\345\255\230\350\260\203\346\225\264", 0, QApplication::UnicodeUTF8));
        exportMeshButton->setText(QApplication::translate("Form", "\345\257\274\345\207\272mesh", 0, QApplication::UnicodeUTF8));
        importSetButton_2->setText(QApplication::translate("Form", "\345\257\274\345\205\245\345\261\200\351\203\250\345\235\220\346\240\207", 0, QApplication::UnicodeUTF8));
        defomationTemButton->setText(QApplication::translate("Form", "\346\250\241\346\235\277\345\275\242\345\217\230", 0, QApplication::UnicodeUTF8));
        repairButton->setText(QApplication::translate("Form", "\346\250\241\346\235\277\347\274\235\345\220\210", 0, QApplication::UnicodeUTF8));
        deformationButton->setText(QApplication::translate("Form", "\345\275\242\345\217\230(\346\212\225\345\275\261)", 0, QApplication::UnicodeUTF8));
        deformationAutoButton->setText(QApplication::translate("Form", "\345\275\242\345\217\230\357\274\210\344\270\200\346\255\245\345\210\260\344\275\215\357\274\211", 0, QApplication::UnicodeUTF8));
        importSetButton->setText(QApplication::translate("Form", "\345\257\274\345\205\245\345\261\200\351\203\250\345\235\220\346\240\207", 0, QApplication::UnicodeUTF8));
        projectBarycentreButton->setText(QApplication::translate("Form", "\346\212\225\345\275\261(\346\214\207\345\220\221\344\270\255\345\277\203)", 0, QApplication::UnicodeUTF8));
        deformationButton2->setText(QApplication::translate("Form", "\345\275\242\345\217\230(\346\212\225\345\275\261)2", 0, QApplication::UnicodeUTF8));
        getBoundaryTemButton->setText(QApplication::translate("Form", "\346\211\276\345\210\260\350\276\271\347\225\214\357\274\210\346\250\241\346\235\277\357\274\211", 0, QApplication::UnicodeUTF8));
        checkBoxShowPoints->setText(QApplication::translate("Form", "\347\234\213\347\202\271", 0, QApplication::UnicodeUTF8));
        deformationAutoButton2->setText(QApplication::translate("Form", "\345\275\242\345\217\2302\357\274\210\344\270\200\346\255\245\345\210\260\344\275\215\357\274\211", 0, QApplication::UnicodeUTF8));
        smoothRootButton->setText(QApplication::translate("Form", "\345\271\263\346\273\221\347\211\231\346\240\271", 0, QApplication::UnicodeUTF8));
        deleteProButton->setText(QApplication::translate("Form", "\345\210\240\351\231\244\346\212\225\345\275\261\351\235\242\347\211\207", 0, QApplication::UnicodeUTF8));
        copyMeshButton->setText(QApplication::translate("Form", "\345\244\215\345\210\266", 0, QApplication::UnicodeUTF8));
        findBoundaryButton->setText(QApplication::translate("Form", "\350\257\206\345\210\253\346\264\236", 0, QApplication::UnicodeUTF8));
        sortBorderButton->setText(QApplication::translate("Form", "\350\276\271\347\225\214\346\216\222\345\272\217", 0, QApplication::UnicodeUTF8));
        fillGapButton->setText(QApplication::translate("Form", "\347\274\235\345\220\210", 0, QApplication::UnicodeUTF8));
        smoothGapButton->setText(QApplication::translate("Form", "\345\271\263\346\273\221\346\216\245\347\274\235", 0, QApplication::UnicodeUTF8));
        cutButton->setText(QApplication::translate("Form", "\345\210\207\345\211\262", 0, QApplication::UnicodeUTF8));
        mergeButton->setText(QApplication::translate("Form", "\345\220\210\345\271\266\346\250\241\345\236\213", 0, QApplication::UnicodeUTF8));
        findNeastestPointsButton->setText(QApplication::translate("Form", "\346\211\276\346\234\200\350\277\221\347\202\271", 0, QApplication::UnicodeUTF8));
        normalizeButton->setText(QApplication::translate("Form", "\346\250\241\345\236\213\346\240\207\345\207\206\345\214\226", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Form", "Step", 0, QApplication::UnicodeUTF8));
        checkBoxHideAxis->setText(QApplication::translate("Form", "\351\232\220\350\227\217\345\261\200\351\203\250\345\235\220\346\240\207", 0, QApplication::UnicodeUTF8));
        exportPButton->setText(QApplication::translate("Form", "\345\257\274\345\207\272\347\202\271", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Form", "Initial", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("Form", "Cut Times", 0, QApplication::UnicodeUTF8));
        toothinfo->setText(QString());
        exportPButton2->setText(QApplication::translate("Form", "\345\257\274\345\207\272\347\202\271\351\233\206", 0, QApplication::UnicodeUTF8));
        errorMapButton->setText(QApplication::translate("Form", "\350\257\257\345\267\256\345\210\206\346\236\220", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("Form", "else Red", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Form", "Blue <=", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("Form", "Cyan <=", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("Form", "Green <=", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("Form", "Yellow <=", 0, QApplication::UnicodeUTF8));
        errorinfo->setText(QString());
        errorinfo_2->setText(QString());
        invertNormalButton->setText(QApplication::translate("Form", "\347\277\273\350\275\254\347\202\271\344\272\221\346\263\225\345\220\221", 0, QApplication::UnicodeUTF8));
        invertNormalFaceButton->setText(QApplication::translate("Form", "\347\277\273\350\275\254mesh\346\263\225\345\220\221", 0, QApplication::UnicodeUTF8));
        extendBorderButton->setText(QApplication::translate("Form", "\346\211\251\345\261\225\350\276\271\347\225\214", 0, QApplication::UnicodeUTF8));
        linkPointsButton->setText(QApplication::translate("Form", "\350\277\236\347\202\271\350\241\245\347\274\235", 0, QApplication::UnicodeUTF8));
        importInsectionPButton->setText(QApplication::translate("Form", "\345\257\274\345\205\245\344\272\244\347\202\271", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Form: public Ui_Form {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_REPAIRDIALOG_H
