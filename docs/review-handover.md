# レビュー引き継ぎ (2026-09-21)

このファイルは一時的なもので、レビューが済んだら消してよい。設計上の判断は
[`two-field-perturbation-design.md`](two-field-perturbation-design.md)、作業の
経緯は [`../daily/2026-09-21.md`](../daily/2026-09-21.md) に書いてある。ここに
あるのは**どの順で何を見るか**だけ。

## 何が積まれているか

`main` (`571b650`) の上に14本、17コミット。二場の摂動スペクトルを、複素modeの
導入から図の生成まで通し、そのあと横断レビューで出た不具合を直したもの。

各ブランチは単独で `cargo test` が緑、`cargo fmt --check` と
`cargo clippy --all-targets` も通ることを確認済み。fast-forwardで下から順に
積める。

## 積む順序

### 機能 (設計docのroadmap 7-12)

| ブランチ | tip | 中身 |
|---|---|---|
| `feature/complex-mode-values` | `4413644` | mode値を `Quantity<D, Complex64>` へ。`normalized_mode_error` |
| `feature/power-spectra` | `a9421b9` | `spectra` module。covariance、射影、reality/range検査 |
| `feature/adaptive-system` | `428fc91` | step制御を `integrate` へ切り出し、4つのsystemが実装 |
| `feature/spectrum-run` | `5e3e5af` | 背景survey、pivot、k gridのdriver、CSV。共形時間の方針をdocに記録 |
| `feature/two-field-spectrum-plots` | `1d38299` | 四パネル図と `generate_all.py` への組み込み |
| `feature/rescaled-mode-variable` | `bdc052e` | `v = a U` を積分する形へ。開始点をe-fold指定に |
| `feature/per-component-tolerances` | `698bb5c` | 絶対許容誤差を成分群ごと・次元付きに |
| `feature/spectrum-validation` | `3ce2098` | 設計docの検証項目のうち未実装やった4つ |

### レビュー由来の修正

| ブランチ | tip | 中身 |
|---|---|---|
| `fix/output-before-horizon-exit` | `8b79ad9` | **地平線を出てへんmodeの真空をスペクトルとして返しとった** |
| `feature/curved-field-space` | `5dcba90` | 曲がった場空間が一度も走ってへんかった |
| `fix/dimensioned-component-error` | `85ec285` | 実数側の誤差関数が素の `f64` を5つ取っとった |
| `fix/fifth-order-crossing` | `aebabf1` | 交差の精密化が誤差制御を通らん4次やった |
| `fix/mode-error-normalization` | `b9096e1` | mode誤差の正規化が根拠と食い違うてた |
| `fix/cross-spectrum-imaginary-part` | `39fb835` | `Im P^12` を無検査で捨てとった |

`docs/review-handover` はこのファイルだけを載せた一番上のブランチ。中身を
読み終えたら落としてよい。

## 特に見てほしい箇所

**`src/quantity.rs` の `PartialEq`/`PartialOrd`。** deriveから手書きに変えた。
deriveは次元markerにまで比較を要求するが、markerはzero-sizedで、そもそも同じ
次元同士しか比較できひん。crateのcore型に触っとるのでここが一番レビューの
価値がある。なお従来のderive版は `MassPower` が `PartialOrd` をderiveして
へんため**実体化不可能**やった (`Time < Time` が書けへん状態)。

**`src/perturbation.rs` の `v = a U` への変換。** 新旧の定式化でspectraが
10桁一致することを確認してある。drag項、`a''/a` 項、Wronskianの前因子、
BD初期化が全部違うのに同じ答えになるので相互検証としては強いが、式そのものは
見てほしい。

**`src/spectrum_run.rs` の `advance_until` と `refine_crossing`。** 全段を
イベントとして書き、最後の1stepをbisectionで詰める。イベントの正側で止める
ので、BD初期化に渡す状態が要求深度より浅くならへん。

## まだ直してへん指摘

レビューで挙がったが未着手のもの。重い順。

1. **`ModeAmplitude` のdocが今も「`U^I_A`」。** 実体は `v^I_A`。書き漏れ。
   併せて `mode_covariance(&modes, ScaleFactor::new(1.0))` と書くと `a^2` 倍
   (実行例では `e^110` 程度) ずれるのに型でもruntimeでも止まらん。公開の
   入口を状態を取る形にして、2引数版を `pub(crate)` に落とすのが素直。
2. **`validate()` の穴。** `output_e_folds_before_end` と `initial_dt` を
   検査してへん。NaNや0で `max_steps` 回まわしてから失敗する。
3. **`SpectrumRunError` に `Display`/`Error` がない。** `?` で運べず、
   `MinimumStepReached` がどの段階・どのkか言わへん。
4. **`Dual`/`Jet` に `exp`/`ln`/`cos`/`tanh` がない。** 多項式と平方根だけ
   なので、Starobinsky・natural・alpha-attractorが**書けへん**。
5. **多場側に `slow_roll_initial_state` がない。** 単一場にはある。二場の
   利用者はexampleからコピペする羽目になる。
6. **`e_s` の向きの規約がAPIにもCSVにも書いてない。** `C_RS` の符号は
   その向き次第。手元のモデルでは採点間隔1〜1000でビット一致やったが、
   規約が書かれてへんこと自体が問題。
7. **`spectrum_run.rs` が940行。** driver、イベント探索、射影、CSVが同居。
   `advance_until` は `integrate` にあるべきで、`background.rs` は別方式の
   イベント処理を再実装しとる。
8. **不要な2。** `TwoFieldModes`、`wronskian`、`mode_covariance` などは
   `U^I_A` がN でも正方なのでconst-genericで書けた。`KinematicBasis2`、
   `project`、`r_C` は本物の二場物理なので隔離は正しい。

設計docの「Numerical and physical validation」節の末尾に、既知の限界として
dynamic rangeとBD近似の2点を書いてある。roadmapの13/14はそこから来ている。

## 確認のしかた

```sh
cargo test                      # compile-failを含む全テスト
cargo clippy --all-targets      # 警告ゼロを維持しとる
cargo fmt --check
uv run plots/generate_all.py    # CSVと図を10ファイル生成して検証
```

`cargo test` は4-5秒、`generate_all.py` は11秒ほど。
