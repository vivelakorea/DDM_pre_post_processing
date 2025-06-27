# Abaqus & ParaDiS VTK 변환 통합 파이프라인

본 스크립트(`DDM2para.py`)는 다음과 같은 시뮬레이션 결과를 VTK 형식으로 자동 변환하고, 결과를 압축 파일로 정리해주는 통합 도구입니다.

- **Abaqus `.odb` 파일 → `.vtu` 파일**
- **ParaDiS `restart` 파일 → `.vtp` 파일**
- **모든 결과를 하나의 `.zip` 파일로 압축**

---

## 🔧 실행 환경

- Python 3.7 이상
- **Abaqus** 설치 및 명령줄에서 `abaqus python` 명령 사용 가능해야 함
- [pyevtk](https://github.com/paulo-herrera/PyEVTK) 설치 필요 (ParaDiS 변환용):
  ```bash
  pip install pyevtk
  ```

---

## 📁 폴더 구조 및 스크립트 위치

아래와 같은 폴더 구조에서 이 스크립트를 실행하세요:

```
├── DDM2para.py
├── ABAQUS/
│   └── Job-compression_biXtal_ori_holding.odb
├── restart/
    ├── rs0001.data
    ├── rs0002.data
    ├── ...
```

### ✅ 스크립트(`DDM2para.py`) 위치
- 반드시 **최상위 폴더(root)** 에 위치해야 합니다.
- 즉, `ABAQUS/`, `restart/`와 같은 폴더와 **같은 경로**에 두어야 하며, 그 하위 폴더가 아닙니다.

### ✅ 실행 위치
`DDM2para.py`가 위치한 폴더에서 터미널 또는 명령 프롬프트를 열고 다음과 같이 실행합니다:

```bash
python DDM2para.py full_run --odbFile Job-compression_biXtal_ori_holding.odb --n-frames 10
```

이렇게 하면 `.odb` 파일은 `ABAQUS/`, `.data` 파일들은 `restart/` 폴더에서 자동으로 불러옵니다.


## 🚀 실행 방법

```bash
python DDM2para.py <모드> [옵션]
```

사용 가능한 실행 모드:

---

### 1. `full_run`  
**Abaqus → ParaDiS → 압축(zip)** 까지 전체 파이프라인 수행

```bash
python DDM2para.py full_run --odbFile 모델.odb [--n-frames N] [--no-parallel]
```

| 옵션명            | 설명 |
|-------------------|------|
| `--odbFile`       | 필수. `ABAQUS/` 폴더 내의 `.odb` 파일 이름 |
| `--n-frames`      | 선택. 변환할 프레임 수 (Abaqus와 ParaDiS 둘 다에 적용) |
| `--no-parallel`   | 선택. 병렬 처리를 비활성화하고 직렬로 수행 |

---

### 2. `paradis_only`  
**ParaDiS → VTK 변환만 수행**

```bash
python DDM2para.py paradis_only --paradis-frames 50 [--n-frames N] [--no-parallel]
```

| 옵션명               | 설명 |
|----------------------|------|
| `--paradis-frames`   | 필수. 변환할 ParaDiS restart 프레임 수 |
| `--n-frames`         | 선택. 위 값을 덮어쓰는 프레임 수 (full_run과 동일하게 맞출 때 유용) |
| `--no-parallel`      | 선택. 병렬 처리 비활성화 |

---

## 📦 출력 결과 구조

| 폴더 | 설명 |
|------|------|
| `ABAQUS/모델이름/` | `.vtu` 파일들과 `.pvd` 파일 저장 |
| `vtk_paradis/`     | `.vtp` 파일들과 `paradis_series.pvd` 저장 |
| `zip_archives/`    | `.vtu`, `.vtp`, `.pvd` 파일들을 묶은 압축파일 저장 |

---

## 🎯 옵션 사용 예시

### 예시 1 – 전체 파이프라인 실행 (프레임 20개만 등간격 추출)
```bash
python DDM2para.py full_run --odbFile mymodel.odb --n-frames 20
```

### 예시 2 – ParaDiS 결과만 100개 프레임 추출 (병렬 비활성화)
```bash
python DDM2para.py paradis_only --paradis-frames 100 --no-parallel
```

---

## 📝 추가 참고 사항

- 프레임 선택은 `np.linspace(start, end, N)` 방식으로 **균일 간격 등분 추출**됩니다.
- `.pvd` 파일이 생성되므로 ParaView 등에서 시간 기반 시각화에 사용 가능합니다.
- 실행 중 `temp_abaqus_worker.py` 라는 임시 스크립트가 생성되며, 실행 종료 후 자동 삭제됩니다.
- `.odb`, `rs****.data` 파일 등은 각각 `ABAQUS/`, `restart/` 폴더에 있어야 합니다.

---

## 📬 문의 및 기여

본 도구는 연구용 스크립트입니다. 내부 구조 개선, 기능 추가, 에러 리포트 등 언제든 환영합니다.

---

