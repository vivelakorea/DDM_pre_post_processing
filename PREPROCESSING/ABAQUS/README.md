# Final INP Generator for Abaqus

이 스크립트는 **ABAQUS CAE에서 생성한 `.inp` 파일을 후처리**하여, 중심 정렬된 모델과 표면 요소(`surface element`), 경계 노드셋(`nodeset`) 정의를 포함한 최종 `.inp` 파일을 생성합니다. 또한 VUMAT 또는 ELA 테스트용 시뮬레이션을 위한 footer를 자동으로 붙여줍니다.

> ⚠️ **입력 파일은 반드시 ABAQUS CAE에서 다음 절차에 따라 생성된 `.inp` 파일이어야 합니다:**
> 1. 파트(Part)를 정의하고,
> 2. 해당 파트를 메싱(Meshing)한 뒤,
> 3. `.inp` 파일로 `Write Input` 을 통해 저장
>
> → 수동으로 작성하거나, 중간 단계에서 임의로 수정된 `.inp` 파일은 예상치 못한 동작을 유발할 수 있습니다.

---

## 🔧 기능 요약

- 모델의 노드 좌표를 중심 정렬 (centered alignment)
- 표면 요소(`*Elset, elset=Set-SurfaceEl`) 자동 추출
- 모델 경계면의 노드셋(`*Nset`) 자동 생성 (X0, XL, Y0, YL, Z0, ZL 등)
- 전체 요소 셋(`*Elset, elset=elset_all`)과 섹션(`*Solid Section`) 추가
- 테스트 조건(예: ELA 또는 VUMAT)에 따라 footer 파일을 자동 병합

---

## 📁 폴더 구조

```
project_root/
│
├── input/
│   ├── input_file.inp         ← 원본 Abaqus 모델 파일
│   └── footer/                ← 필요한 footer 조각 파일
│       ├── footer1_VUMAT.inp
│       ├── footer2_VUMAT.inp
│       └── ...
│
├── output/
│   └── Job-final_output.inp   ← 생성된 최종 파일 (자동 저장됨)
│
└── generate_final_inp.py      ← 이 Python 스크립트 파일
```

---

## 🖥️ 사용법

```bash
python generate_final_inp.py --input your_model.inp --output final_output.inp
```

옵션 설명:

- `--input`: `input/` 폴더 안의 입력 `.inp` 파일 이름 (확장자 `.inp`는 생략해도 됨)
- `--output`: 출력 `.inp` 파일 이름 (자동으로 `output/Job-<이름>.inp`에 저장됨)
- `--ela_test`: ELA 테스트용 footer 사용 시 이 플래그 추가

예시:

```bash
python generate_final_inp.py --input beam_model --output beam_ready
```

ELA 테스트용:

```bash
python generate_final_inp.py --input beam_model --output beam_ela --ela_test
```

---

## 📌 주의사항

- 입력 파일은 반드시 `input/` 폴더에 있어야 합니다.
- footer 파일들은 `input/footer/` 폴더에 있어야 하며, 다음 중 하나의 세트를 구성해야 합니다:
  - VUMAT용: `footer*_VUMAT.inp`
  - ELA용: `footer*_ELA.inp`
- 출력은 자동으로 `output/Job-*.inp`에 저장됩니다.
- Abaqus 모델 내 파트 이름은 `"Part-1"`으로 되어 있어야 하며, 자동으로 `"Microbeam"`으로 대체됩니다.

---

## 🧱 생성되는 Nset 이름 규칙

모서리 경계면 노드셋은 다음과 같은 규칙으로 자동 생성됩니다:

- `X0`: 최소 x값에 있는 노드들
- `XL`: 최대 x값에 있는 노드들
- `Y0`, `YL`, `Z0`, `ZL`: 각각 y, z 방향에 대해 동일

복합 경계면도 자동 생성됩니다:

- 예: `X0Y0ZL` → x=최소, y=최소, z=최대인 노드들

---

## 📞 문의

작성자: Gyujang Sim  
소속: MAMEL Lab, SNU  
문의: gyujang95@snu.ac.kr
