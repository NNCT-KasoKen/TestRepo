//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Geometry.cc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "Geometry.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleRunge.hh"
#include "G4AutoDelete.hh"

//磁場系のやつ
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4VUserDetectorConstruction.hh"

//電場系の奴
#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"


//ジオメトリを引き算するためのhh
#include "G4SubtractionSolid.hh"







/*
//自分で磁場頑張るやつ
#include "G4PhysicalConstants.hh"
class MyMagneticField : public G4MagneticField
{
public:
    MyMagneticField() : G4MagneticField(), fB0(0 * tesla), fGradient(1000 / (1.0 * second)) {}

    virtual void GetFieldValue(const G4double point[4], G4double *Bfield) const
    {
        // point[0], point[1], point[2] は位置座標(x, y, z)で、point[3] は時間(t)です。
        G4double time = point[3];

        // 1次関数的に変化する磁場を設定
        G4double Bz = fB0 + fGradient * time;

        // 磁場はz方向にのみ影響を与えるものと仮定
        Bfield[0] = 0.0;
        Bfield[1] = Bz;
        //Bfield[1] = 5* tesla;
        Bfield[2] = 0.0;
    }

private:
    G4double fB0;        // 初期磁場の強さ
    G4double fGradient;  // 時間に関する勾配
};
*/







Geometry::Geometry() {}

Geometry::~Geometry() {}

G4VPhysicalVolume* Geometry::Construct(){
  G4NistManager* materi_Man = G4NistManager::Instance();

  //世界の大きさと形を決める
  G4double leng_X_World = 5.0*m;
  G4double leng_Y_World = 5.0*m;
  G4double leng_Z_World = 5.0*m;
  auto solid_world =new G4Box{"Solid_WOrld", leng_X_World/2.0, leng_Y_World/2.0, leng_Z_World/2.0};//値まとめ


  //---------空気作り---------
  //世界の材質を設定する
  //まずは自分でいい感じな気圧の空気を作る
  //~~原子作り~~
  //水素原子
  G4double hydrogen_a = 1.008*g/mole;
  G4Element* elH = new G4Element("Hydrogen", "H",1.,hydrogen_a);
  //ヘリウム原子
  G4double helium_a = 4.003*g/mole;
  G4Element* elHe = new G4Element("Helium", "He",2.,helium_a);
  //炭素原子
  G4double carbon_a = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon", "C",6.,carbon_a);
  //窒素原子
  G4double nitrogen_a = 14.007*g/mole;
  G4Element* elN = new G4Element("Nitrogen", "N",7.,nitrogen_a);
  //酸素原子
  G4double oxygen_a = 16.000*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",8.,oxygen_a);
  //ネオン原子
  G4double neon_a=20.180*g/mole;
  G4Element* elNe = new G4Element("Neon","Ne",10.,neon_a);
  //アルゴン原子
  G4double argon_a=39.948*g/mole;
  G4Element* elAr = new G4Element("Argon","Ar",18.,argon_a);
  //クルプトン原子
  G4double krypton_a=83.80*g/mole;
  G4Element* elKr = new G4Element("Krypton","Kr",36.,krypton_a);

  //~~化合物作り~~
  //窒素分子
  G4double density_N2 = 1.2507*mg/cm3;
  G4Material* N2 = new G4Material("NitrogenMolecule",density_N2,1,kStateGas);
  N2->AddElement(elN, 2);
  //酸素分子
  G4double density_O2=1.4289*mg/cm3;
  G4Material* O2 = new G4Material("OxygenMolecule",density_O2,1,kStateGas);
  O2->AddElement(elO,2);
  //アルゴン単体
  G4double density_Ar=1.7828*mg/cm3;
  G4Material* Ar = new G4Material("Argon",density_Ar,1,kStateGas);
  Ar->AddElement(elAr,1);
  //二酸化炭素
  G4double density_CO2=1.9768*mg/cm3;
  G4Material* CO2 = new G4Material("CarbonDioxide",density_CO2,2,kStateGas);
  CO2->AddElement(elC,1);
  CO2->AddElement(elO,2);
  //一酸化炭素
  G4double density_CO=1.2501*mg/cm3;
  G4Material* CO = new G4Material("CarbonMonoxide",density_CO,2,kStateGas);
  CO->AddElement(elC,1);
  CO->AddElement(elO,1);
  //ネオン単体
  G4double density_Ne=0.8713*mg/cm3;
  G4Material* Ne = new G4Material("Neon",density_Ne,1,kStateGas);
  Ne->AddElement(elNe,1);
  //ヘリウム単体
  G4double density_He=0.1769*mg/cm3;
  G4Material* He = new G4Material("Helium",density_He,1,kStateGas);
  He->AddElement(elHe,1);
  //メタン(CH4)
  G4double density_CH4=0.717*mg/cm3;
  G4Material* CH4 = new G4Material("Methane",density_CH4,2,kStateGas);
  CH4->AddElement(elC,1);
  CH4->AddElement(elH,4);
  //クルプトン単体
  G4double density_Kr=3.75*mg/cm3;
  G4Material* Kr = new G4Material("Krypton",density_Kr,1,kStateGas);
  Kr->AddElement(elKr,1);
  //一酸化二窒素
  G4double density_N2O=2.8122*mg/cm3;
  G4Material* N2O = new G4Material("DinitrogenOxide",density_N2O,2,kStateGas);
  N2O->AddElement(elN,2);
  N2O->AddElement(elO,1);
  //水素分子
  G4double density_H2 = 0.0898*mg/cm3;
  G4Material* H2 = new G4Material("HydrogenMolecule",density_H2,1,kStateGas);
  H2->AddElement(elH, 2);
  //オゾン
  G4double density_O3=2.141*mg/cm3;
  G4Material* O3 = new G4Material("Ozone",density_O3,1,kStateGas);
  O3->AddElement(elO,3);


  /*既存の空気(G4_AIR)と同じ値にしたかった時のコード*/
  /*
  G4Material* AIR = materi_Man->FindOrBuildMaterial("G4_AIR");
  G4double density_AIR = AIR->GetDensity();
  G4double temperature_AIR = AIR->GetTemperature();
  G4double pressure_AIR = AIR->GetPressure();
  */

  //空気の数値の定義　PV=nRTより、密度はn/V=P/TR
  G4double pressure_AIR = 1*pascal;
  G4double temperature_AIR = 300*kelvin;

  //密度を計算で出す。
  G4double gas_constant = 8.31;
  G4double density_AIR = 2*(pressure_AIR / (temperature_AIR*gas_constant));
  //G4double density_AIR = 1562.5*g/m3;
  
  //実際に1つの混合物(空気)にする
  //G4Material* materi_world = new G4Material("My_AIR", density_AIR, 12,kStateGas,temperature_AIR,pressure_AIR);
  G4Material* materi_world = new G4Material("My_AIR", density_AIR, 1,kStateGas,temperature_AIR,pressure_AIR);
  //G4Material* materi_world = materi_Man->FindOrBuildMaterial("G4_Galactic");//真空
  //G4Material* materi_world = materi_Man->FindOrBuildMaterial("G4_AIR");//空気

  /*容積比
  materi_world->AddMaterial(N2,78.11*perCent);
  materi_world->AddMaterial(O2,20.96*perCent);
  materi_world->AddMaterial(Ar,0.9343*perCent);
  materi_world->AddMaterial(CO2,0.03*perCent);
  materi_world->AddMaterial(CO,1.e-3*perCent);
  materi_world->AddMaterial(Ne,1.8e-3*perCent);
  materi_world->AddMaterial(He,5.3e-4*perCent);
  materi_world->AddMaterial(CH4,1.52e-4*perCent);
  materi_world->AddMaterial(Kr,1.e-4*perCent);
  materi_world->AddMaterial(N2O,5.e-5*perCent);
  materi_world->AddMaterial(H2,5.e-5*perCent);
  materi_world->AddMaterial(O3,2.e-5*perCent);
  */
  /*重量比*/
  /*
  materi_world->AddMaterial(N2,75.53*perCent);
  materi_world->AddMaterial(O2,23.14*perCent);
  materi_world->AddMaterial(Ar,1.280*perCent);
  materi_world->AddMaterial(CO2,0.045*perCent);
  materi_world->AddMaterial(CO,1.e-5*perCent);
  materi_world->AddMaterial(Ne,1.2e-5*perCent);
  materi_world->AddMaterial(He,7.3e-5*perCent);
  materi_world->AddMaterial(CH4,8.4e-3*perCent);
  materi_world->AddMaterial(Kr,3.e-4*perCent);
  materi_world->AddMaterial(N2O,8.e-5*perCent);
  materi_world->AddMaterial(H2,3.e-6*perCent);
  materi_world->AddMaterial(O3,3.e-6*perCent);
  */
  materi_world->AddMaterial(H2,100*perCent);

  auto logVol_world = new G4LogicalVolume{solid_world, materi_world, "LogVol_World"};//まとめ
  logVol_world->SetVisAttributes(G4VisAttributes::GetInvisible());//視覚情報の設定、世界は見えなくする


  //世界を置く
  G4int copyNum_world = 0;//ID的な何か
  auto physVol_world = new G4PVPlacement{G4Transform3D(), "PhysVol_World", logVol_world, 0, false, copyNum_world};//実際に世界を作る(今までの設定をまとめる)




  //---------ジオメトリ---------
  //ジオメトリの材質を決める
  G4Material* geoMateri = materi_Man->FindOrBuildMaterial("G4_Fe");//鉄
  //G4Material* geoMateri = materi_Man->FindOrBuildMaterial("G4_BRASS");//真鍮(実際に使われる素材)

  //物体のを作る(大小2角円柱を引き算して作る)
  //大きいほうの円柱
  G4double innerRadiusBig = 0*cm;//内側の半径
  G4double outerRadiusBig = 1.9*cm;//外側の半径 空間の半径15mm+厚さ0.4mm
  G4double lengZBig = 2.8*cm;//円柱の高さ(厚さ) 上の厚さ4mm+空間の厚さ20mm+下の厚さ4mm=24mm
  G4double startDegreeBig = 0.*deg;//始まりの角度
  G4double endDegreeBig = 180.*deg;//終端角度
  auto tubeBig = new G4Tubs{"TUBE_BIG", innerRadiusBig*2, outerRadiusBig*2, lengZBig, startDegreeBig, endDegreeBig};

  //小さいほうの円柱
  G4double innerRadiusSmall = 0*cm;//内側の半径
  G4double outerRadiusSmall = 1.5*cm;//外側の半径
  G4double lengZSmall = 2.0*cm;//円柱の高さ(厚さ)
  G4double startDegreeSmall = 0.*deg;//始まりの角度
  G4double endDegreeSmall = 180.*deg;//終端角度
  auto tubeSmall = new G4Tubs{"TUBE_SMALL", innerRadiusSmall*2, outerRadiusSmall*2, lengZSmall, startDegreeSmall, endDegreeSmall};

  //ジオメトリの引き算操作
  auto cmpGeo = new G4SubtractionSolid{"cmpGeo",tubeBig,tubeSmall};
  
  //大きさ(形)と材質をまとめる。物体の完成
  auto logVol = new G4LogicalVolume{cmpGeo, geoMateri, "LogVol_cmpGeo", 0, 0, 0};



  //視認性向上のためにジオメトリに色を付与
  G4VisAttributes* color = new G4VisAttributes(G4Colour(1, 0, 0));
  logVol->SetVisAttributes(color);

  

  //位置
  G4double pos_X_LogV = 0.0*cm;
  G4double pos_Y_LogV = 0.0*cm;
  G4double pos_Z_LogV = 0.0*cm;
  G4double gap = 0.30*cm;
  auto threeVect1 = G4ThreeVector{pos_X_LogV, pos_Y_LogV, pos_Z_LogV+(gap/2)};
  auto threeVect2 = G4ThreeVector{pos_X_LogV, pos_Y_LogV, pos_Z_LogV-(gap/2)};
  //角度
  auto rotMtrx1 = G4RotationMatrix();
  auto rotMtrx2 = G4RotationMatrix();
  rotMtrx1.rotateX(90.*deg);
  rotMtrx2.rotateX(90.*deg);
  rotMtrx2.rotateY(180.*deg);

  //位置と角度まとめ
  auto trans3D1 = G4Transform3D{rotMtrx1, threeVect1};
  auto trans3D2 = G4Transform3D{rotMtrx2, threeVect2};

  //ID
  G4int copyNum1 = 10;
  G4int copyNum2 = 20;

  //実際に置く
  new G4PVPlacement{trans3D1, "PhysVol_TUBE1", logVol, physVol_world, false, copyNum1, true};
  new G4PVPlacement{trans3D2, "PhysVol_TUBE2", logVol, physVol_world, false, copyNum2, true};



  //---------磁場かけるためのための実体のないジオメトリ---------a
  //ジオメトリの材質を決める
  //G4Material* geoMagMateri = materi_Man->FindOrBuildMaterial("G4_Galactic");//実体をなくすので真空
  //G4Material* geoMagMateri = materi_Man->FindOrBuildMaterial("G4_AIR");//世界と同じがいいかもね


  //logVolの設定
  //auto logVol_Mag = new G4LogicalVolume{tubeBig, geoMagMateri, "LogVol_Mag", 0, 0, 0};
  auto logVol_Mag = new G4LogicalVolume{tubeBig, materi_world, "LogVol_Mag", 0, 0, 0};

  //磁場の設定　既存の磁場を使っている。既存のを使うと直線的な磁場になってしまう
  G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(0.,5.0*tesla,0.));
  G4FieldManager* localFieldMgr_mag = new G4FieldManager(magField);
  logVol_Mag->SetFieldManager(localFieldMgr_mag,true);//磁場を付与

  //視認性向上のためにジオメトリに色を付与(真空だけど色がつく)。G4VisAttributesの第1引数が可視属性、第2で色を設定できる。
  G4VisAttributes* color_Mag = new G4VisAttributes(false,G4Colour(0, 0, 0));
  logVol_Mag->SetVisAttributes(color_Mag);
  

  //ID
  G4int copyNumMag1 = 110;
  G4int copyNumMag2 = 120;

  //実際に置く
  new G4PVPlacement{trans3D1, "PhysVol_Mag1", logVol_Mag, physVol_world, false, copyNumMag1, true};
  new G4PVPlacement{trans3D2, "PhysVol_Mag2", logVol_Mag, physVol_world, false, copyNumMag2, true};


  //---------電場かけるためのための実体のないジオメトリ---------
  //ジオメトリの大きさと形を決める
  //ジオメトリが重なってると処理がおかしくなるかもしれない
  //→x方向の大きさとギャップを微妙に小さくしている。
  G4double leng_X_EleField = outerRadiusBig*2 -1.*nm;
  G4double leng_Y_EleField = lengZBig*2;
  G4double leng_Z_EleField = gap-1.*nm;
  auto solid_EleField =new G4Box{"Solid_WOrld", leng_X_EleField/2.0, leng_Y_EleField/2.0, leng_Z_EleField/2.0};//値まとめ

  //logVolにまとめる
  auto logVol_EleField1 = new G4LogicalVolume{solid_EleField, materi_world, "LogVol_EleField", 0, 0, 0};
  auto logVol_EleField2 = new G4LogicalVolume{solid_EleField, materi_world, "LogVol_EleField", 0, 0, 0};


  //電場の設定
  
  G4ElectricField* eleField = new G4UniformElectricField(G4ThreeVector(0.5 * volt / cm, 0.,0. ));
  G4FieldManager* localFieldMgr_ele = new G4FieldManager(eleField);
//  logVol_EleField1->SetFieldManager(localFieldMgr_ele, true);  // 電場を付与
//  logVol_EleField2->SetFieldManager(localFieldMgr_ele, true);
  

  //視認性向上のためにジオメトリに色を付与(真空だけど色がつく)。G4VisAttributesの第1引数が可視属性、第2で色を設定できる。
//  logVol_EleField1->SetVisAttributes(new G4VisAttributes(true,G4Colour(1, 0.573, 0.286)));
//  logVol_EleField2->SetVisAttributes(new G4VisAttributes(true,G4Colour(0.796, 0.251, 0.871)));

  //ID
  G4int copyNumEleField = 210;

  //位置
  G4double pos_X_EleField = outerRadiusBig;
  G4double pos_Y_EleField = 0.0*cm;
  G4double pos_Z_EleField = 0.0*cm;
  
  auto threeVectEleField = G4ThreeVector{pos_X_EleField, pos_Y_EleField, pos_Z_EleField};
  //角度
  auto rotMtrxEleField = G4RotationMatrix();

  //位置と角度まとめ
  auto trans3DEleField = G4Transform3D{rotMtrxEleField, threeVectEleField};

  //実際に置く
  //1つ目
  new G4PVPlacement{trans3DEleField, "PhysVol_EleField", logVol_EleField1, physVol_world, false, copyNumEleField, true};
  //2つ目
  threeVectEleField = G4ThreeVector{pos_X_EleField*-1, pos_Y_EleField, pos_Z_EleField};
  rotMtrxEleField.rotateX(180.*deg);
  trans3DEleField = G4Transform3D{rotMtrxEleField, threeVectEleField};
  copyNumEleField = 220;
  new G4PVPlacement{trans3DEleField, "PhysVol_EleField", logVol_EleField2, physVol_world, false, copyNumEleField, true};

  //世界を返す
  return physVol_world;
}
