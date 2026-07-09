# sentinel-findings.md

*Written by Sentinel — WARN-tier findings.*

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:33.432419+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1099`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:33.432460+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1101`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:33.432475+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1104`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:33.432486+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1106`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:11:33.432496+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:11.261461+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1099`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:11.261488+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1101`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:11.261497+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1104`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:11.261504+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1106`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:33:11.261510+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:19.824259+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1099`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:19.824293+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1101`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:19.824306+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1104`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:19.824315+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1106`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:04:19.824324+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1098`
- **Detail:** pattern matched: analysis_name = str(clean["Analysis.name"].iloc[0]) if "Analysis.name" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:26:01.559055+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1099`
- **Detail:** pattern matched: subgroup = str(clean["Subgroup"].iloc[0]) if "Subgroup" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:26:01.559089+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1101`
- **Detail:** pattern matched: str(clean["Analysis.group"].iloc[0]) if "Analysis.group" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:26:01.559101+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1104`
- **Detail:** pattern matched: str(clean["Analysis.number"].iloc[0]) if "Analysis.number" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:26:01.559110+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `scripts/run_pairwise70_benchmark.py:1106`
- **Detail:** pattern matched: review_doi = str(clean["review_doi"].iloc[0]) if "review_doi" in clean.columns else ""
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T13:26:01.559118+00:00

## [WARN] P1-license-noncompliance
- **Location:** `(repo-wide)`
- **Detail:** package 'pyreadr' has license 'GNU Affero General Public License v3 or later (AGPLv3+)' which conflicts with MIT-licensed publication. If this repo is intentionally GPL/AGPL, suppress with a sentinel:skip-line marker; otherwise swap to a permissive-licensed alternative.
- **Fix hint:** Options: (1) replace the dep with an MIT/BSD/Apache-2.0 alternative; (2) change the repo's LICENSE to GPL/AGPL to match the dep; (3) if dual-licensed in your favor, add a comment near the dep declaration with the chosen license and add 'sentinel:skip-line P1-license-noncompliance' nearby.
- **Source:** lessons.md#license-compliance-2026-04-29
- **When:** 2026-05-15T15:39:37.124235+00:00
