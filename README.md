# DDD-Bicrystal Toolkit Repository

본 저장소는 ParaDiS와 Abaqus를 기반으로 하는 **DDM 시뮬레이션과 후처리 파이프라인**을 다룹니다. 아래는 각각의 서브 모듈 및 스크립트에 대한 설명입니다.

---

## 📁 구성 요소 개요

### 1. `ParaDiS Bicrystal Generator (MATLAB)`
- **위치**: `bicrystal_generator/`
- **기능**: Frank–Read 소스를 기반으로 bicrystal 초기 조건 생성 (`.data`, `.ctrl`)
- **주요 특징**:
  - grain boundary가 `y = 0`에 고정된 bicrystal
  - y > 0 영역의 결정은 Bunge Euler 각으로 회전 가능
  - 목표 전위 밀도 달성 시까지 무작위 소스 배치
  - ParaDiS 실행용 초기 상태 생성

📄 자세한 설명: [`README.md`](./POSTPROCESSING/README.md)

---

### 2. `Final INP Generator (Python)`
- **위치**: `inp_postprocess/`
- **기능**: ABAQUS CAE로부터 생성된 `.inp` 파일 후처리 및 변환
- **주요 특징**:
  - 노드 중심 정렬
  - 표면 요소 자동 추출
  - 경계 노드셋(`X0`, `XL`, `Y0`, `ZL` 등) 자동 생성
  - footer 파일을 붙여 `.inp` 완성

📄 자세한 설명: [`README.md`](./PREPROCESSING/ABAQUS/README.md)

---

### 3. `VTK 변환 파이프라인 (Python)`
- **위치**: `vtk_pipeline/`
- **기능**: ParaDiS 및 Abaqus 결과를 `.vtu`, `.vtp`로 변환 + zip 압축
- **주요 특징**:
  - Abaqus `.odb` → `.vtu`
  - ParaDiS `restart` → `.vtp`
  - `.pvd` 포함 시계열 시각화 가능
  - 전체 자동화 및 병렬 처리 지원

📄 자세한 설명: [`README.md`](./PREPROCESSING/PARADIS/bi_crystal/README.md)

---

## 📞 문의

작성자: Gyujang Sim  
소속: MAMEL Lab, SNU  
문의: gyujang95@snu.ac.kr
