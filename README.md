# 🧬 Full Spectrum Bio-Analysis Pipeline & Monitoring Dashboard

> **Full Spectrum 생의학 데이터 분석 파이프라인 및 실시간 모니터링 대시보드**  
> 유전체 데이터(VCF) 및 임상 데이터의 자동 수집, 정제, 다차원 분석(PRS 계산, 혼합 효과 회귀 분석), 데이터베이스 적재, 그리고 이를 실시간 모니터링할 수 있는 플러머(Plumber) REST API 및 대시보드 웹 서비스가 결합된 통합 생물정보학(Bioinformatics) 플랫폼입니다.

---

## 🛠️ 시스템 아키텍처 및 데이터 흐름

```mermaid
graph TD
    subgraph 1. Data Ingestion & Orchestration (Snakemake)
        VCF[Dataset/Genotypes.vcf.gz] -->|bcftools annotate| AVCF[Annotated.vcf.gz]
        GWAS[GWAS Catalog TSV] -->|Extract weights| W[gwas_weights.txt]
        AVCF & W -->|PLINK 1.9| PRS[prs_out.profile]
        
        RawClin[health_raw.csv] -->|clean_health_data.R| CleanClin[health_cleaned.csv]
    end

    subgraph 2. Modeling & Storage
        PRS & CleanClin -->|run_modeling.R| Models[Linear Mixed Model & lm]
        CleanClin & Models -->|load_to_db.R| DB[(SQLite bio_analysis.db)]
    end

    subgraph 3. Presentation Layer (Port 8085)
        DB -->|R Plumber API| API[app.R REST API]
        API -->|HTML/CSS/JS + Mermaid.js| Dash[System Dashboard]
        API -->|R Markdown| Rep[Automated Patient Report]
    end
    
    classDef pass fill:#2ecc71,stroke:#27ae60,color:#fff;
    classDef process fill:#3498db,stroke:#2980b9,color:#fff;
    class DB,API,Dash,Rep process;
```

---

## ✨ 핵심 기능

1. **지능형 데이터 수집 및 전처리 (Data Ingestion & Cleaning)**
   * VCF 염색체 22번 데이터에서 SNP 변이 정보를 자동으로 식별하고, 이상치(음수 BMI 등)와 결측치를 포함한 임상 모의 데이터를 정제합니다.
   * 성별 표준화, 중앙값 기반의 결측 BMI 및 콜레스테롤 보정(Imputation) 기능을 수행합니다.
2. **고성능 변이 어노테이션 및 PRS 스코어 계산 (Orchestration & Compute)**
   * `bcftools`를 통해 변이 ID가 없는 VCF 파일에 고유 식별자(`%CHROM:%POS:%REF:%ALT`)를 부여합니다.
   * 고속 유전체 분석 도구인 `PLINK 1.9`를 로컬에서 구동하여 각 환자별 Polygenic Risk Score (PRS)를 병렬 연산합니다.
3. **고급 통계 모델링 (Mixed-Effects Modeling)**
   * 일반 다중 회귀 모델뿐만 아니라, R `nlme` 패키지를 사용하여 연령대(Age Decade)를 변동 효과(Random Effect)로 정의한 **선형 혼합 효과 모델(Linear Mixed-Effects Model)**을 분석하여 임상적 신뢰도를 극대화합니다.
4. **최적화된 데이터웨어하우스 구축 (SQLite DB)**
   * 정규화된 테이블(`health_records`, `prs_results`, `model_coefficients` 등) 구조를 생성하고, `user_id` 및 `run_id`에 다중 복합 인덱스를 설정하여 1ms 이하의 초고속 쿼리 속도를 제공합니다.
5. **실시간 모니터링 대시보드 & REST API (Port 8085)**
   * **Mermaid.js 기반 실시간 DAG 모니터**: 파이프라인의 핵심 파일 8개의 존재 여부를 백엔드 하트비트(`/system_status`)로 실시간 수신하여 성공(Green), 누락(Red) 상태를 대시보드 상에 실시간으로 시각화합니다.
   * **무결성 검증 러너**: 디스크 파일 규격, 데이터 적재 상태, 통계 모듈 유효성을 즉시 평가하는 5단계 백엔드 테스트를 UI에서 버튼 클릭으로 실행합니다.
   * **맞춤형 환자 건강성적표**: 특정 환자를 선택하면 R Markdown을 통해 ggplot2 인구 분포 내 위치와 맞춤형 진단 정보 카드를 포함한 Dynamic HTML 리포트를 생성해 브라우저에 임베드합니다.

---

## 📁 프로젝트 구조

```text
Bio-Analysis/
├── api/
│   ├── app.R               # R Plumber REST API 서버
│   └── dashboard.html      # 모니터링 대시보드 UI (Dark Theme)
├── bin/
│   └── plink               # PLINK 1.9 로컬 실행 바이너리 (Mac/HPC 환경 지원)
├── database/
│   └── bio_analysis.db     # SQLite 표준 관계형 데이터베이스
├── Dataset/                # 원자료 및 파이프라인 산출 데이터
│   ├── Annotated.vcf.gz    # ID 매핑된 유전체 데이터
│   ├── health_cleaned.csv  # 이상치 및 결측 보정된 임상 데이터
│   ├── model_results.RData # 학습 완료된 R 모델 통계 객체
│   └── model_summary.txt   # R 통계 모델 적합도 리포트 요약본
├── docker/
│   ├── Dockerfile          # 분석 및 Plumber 서버 실행용 Docker 이미지 빌드 파일
│   ├── docker-compose.yml  # 파이프라인 구동 및 웹 API 자동 실행 오케스트레이터
│   └── info.txt            # 시스템 세부 구현 및 검증 내역 정리
├── reports/
│   └── report_template.Rmd # 개별 환자 맞춤형 R Markdown 레포트 템플릿
├── scripts/
│   ├── clean_health_data.R # 임상 데이터 전처리 스크립트
│   ├── generate_mock_data.R# GWAS 가중치 및 임상 모의 데이터 빌더
│   ├── load_to_db.R        # SQLite DB 스키마 생성 및 데이터 적재 스크립트
│   ├── run_modeling.R      # 다중회귀 및 혼합효과 모델 분석 스크립트
│   └── run_tests.R         # 5단계 무결성 자동 검증 러너
├── Snakefile               # Snakemake 파일 기반 파이프라인 의존성 정의서
├── boot.txt                # 서비스 가동 퀵스타트 안내 가이드
└── SOLID_Principles.md     # 객체 지향 5대 설계 원칙 참고 가이드
```

---

## 🚀 시작하기

### 방법 1. Docker Compose (추천 - 단일 명령 실행)

컨테이너 기술을 통해 복잡한 R 패키지 및 바이오 도구(PLINK, bcftools 등) 설치 없이 즉시 실행할 수 있습니다.

```bash
# docker 디렉토리로 이동하여 컨테이너 빌드 및 가동
docker compose -f docker/docker-compose.yml up --build -d
```
* **동작:** 파이프라인(`Snakefile` 단계 전체)을 최초 1회 순차 실행하여 DB 및 모델 파일 빌드를 끝마친 후, **`8085` 포트**로 웹 API 및 모니터링 대시보드를 즉시 호스팅합니다.

### 방법 2. 로컬 실행 (Local Runtime)

**사전 필수 설치 도구:** R(>=4.3.0), Python 3(및 Snakemake), bcftools, PLINK 1.9

1. **파이프라인 전체 빌드 (Snakemake)**
   ```bash
   snakemake --cores 1
   ```
2. **R Plumber API 서버 가동**
   ```bash
   Rscript api/app.R
   ```
3. **대시보드 접속**
   웹 브라우저를 열고 다음 주소에 접속합니다:
   * **[http://localhost:8085/dashboard](http://localhost:8085/dashboard)**

---

## 📡 REST API 엔드포인트 규격

| HTTP Method | Endpoint | Description | Query Parameters |
| :--- | :--- | :--- | :--- |
| **GET** | `/dashboard` | 실시간 모니터링 대시보드 웹 UI 서빙 | 없음 |
| **GET** | `/system_status` | 서버 시간, 프로젝트 루트, 핵심 파일 유무 및 DB 통계 반환 | 없음 |
| **GET** | `/samples` | 데이터베이스에 등록된 유효 환자 ID 목록 조회 | 없음 |
| **GET** | `/run_tests` | 5단계 시스템 무결성 자동 테스트 원격 실행 및 결과 반환 | 없음 |
| **GET** | `/get_prs` | 환자 1명의 임상 수치, PRS 위험군 정보 및 회귀 모델 예측 반환 | `user_id` (예: `HG00096`) |
| **GET** | `/report` | 환자 1명의 동적 시각화 건강성적표 HTML 파일 생성 및 반환 | `user_id` (예: `HG00096`) |
| **POST** | `/run_pipeline` | 콘솔 백그라운드에서 Snakemake 파이프라인 전체를 강제 재실행 | 없음 |
