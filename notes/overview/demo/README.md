# phys-type-demo

`phys-type-overview.tex` のデモ。式を書き間違えるとコンパイルが通らないことを、実際に `cargo build` して見せる。

リポジトリの workspace には入っていない独立した crate なので、このディレクトリで実行する。

```sh
cd notes/overview/demo
```

初回はビルドに少し時間がかかるので、見せる前に一度 `cargo run` を走らせておく。`cargo build --examples` は壊れた例も含むので失敗する。

## 動く例

| コマンド | 内容 | 出力 |
| --- | --- | --- |
| `cargo run` | Friedmann 方程式 $H=\sqrt{(\dot\phi^2/2+V)/(3M_{\rm pl}^2)}$ | `H = 5.787918e-6 M_pl` |
| `cargo run --example slow_roll` | $\epsilon_H$ を Planck 単位と GeV で | 型つきの式はどちらも `5.000e-3`。$M_{\rm pl}^2$ を落とした手書きの式は GeV で `2.965e34` |
| `cargo run --example jet` | 同じポテンシャル $\lambda\phi^4/4$ を `f64` と `Jet<1>` で評価 | `V = 4`、`V' = 8`、`V'' = 12` |
| `cargo run --example tensor` | 極座標 $G=\mathrm{diag}(1,\rho^2)$ での運動項と力 | $G_{IJ}\dot\phi^I\dot\phi^J=0.37$（手書きは 0.10）、力 $(10^{-3},10^{-3})$（手書きは $\theta$ 成分が 4 倍） |

## コンパイルが通らない例

`cargo build --example <名前>` で、どれもエラーになる。

### 次元

| 名前 | 間違い | 型で見ると |
| --- | --- | --- |
| `forgot_square` | 分母の $M_{\rm pl}^2$ を $M_{\rm pl}$ に | $M^{3/2}$ を `Hubble`（$M$）に代入 |
| `add_velocity` | $\dot\phi^2/2$ の代わりに $\dot\phi$ | $M^2+M^4$ |
| `forgot_planck_mass` | $\epsilon_H$ の $M_{\rm pl}^2$ を落とす | $M^2$ を `Dimensionless` に代入 |
| `conformal_time_step` | 共形時間を宇宙時間として渡す | 次元はどちらも $M^{-1}$。`ConformalTime` は別の型 |
| `pivot_in_mpc` | 観測の $0.05\,\mathrm{Mpc}^{-1}$ を内部の $k$ として渡す | `ComovingWavenumberMpcInverse` は `Mass` ではない |

### テンソル（極座標）

| 名前 | 間違い | 型で見ると |
| --- | --- | --- |
| `forget_metric` | $\dot\phi^I\dot\phi^I$（計量を忘れる） | 上付き同士は `contract` できない |
| `gradient_as_force` | $G^{IJ}V_{,J}$ の代わりに $V_{,I}$ | 上付きと下付きは足せない。次元は合っているので次元の検査だけでは見つからない |
| `raise_with_metric` | $G_{IJ}$ で添字を上げる | 下付き2つの計量には `raise` がない |

## 見せるときの注意

- 次元のエラー文は `typenum` の内部表現（`PInt<UInt<UInt<_, B0>, B0>>` など）が出るので読みにくい。`mismatched types` の行と、例のファイル先頭のコメントを合わせて見せる。
- テンソルのエラー文は `Upper` / `Lower` がそのまま出るので、そのまま見せられる。
- `slow_roll` は「Planck 単位のテストでは気づけない」例。`forgot_planck_mass` と続けて見せる。
- `tensor` は動く例を先に見せ、そのあと `forget_metric`・`gradient_as_force`・`raise_with_metric` を見せる。
