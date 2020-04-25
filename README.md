-------------
practice 1:<br>
  1. 座標変換
    •概要<br>
オイラー角を用いて方向余弦行列を計算し,座標変換を行う. 

practice 2:<br>
  2. Cartesian から keplerian への変換<br>
    •概要<br>
    直交座標系で表現された状態ベクトル(位置+速度)からケプラーの軌道要素を計算する. 

practice 3:<br>
  3. keplerian から Cartesian への変換<br>
    •概要<br>
    ケプラーの軌道6要素から直交座標系における状態ベクトルを計算する.<br>
    具体的には軌道要素から計算した焦点中心座標系での状態ベクトルを座標変換により直交座標系で表現する. 

practice 4:<br>
  4. 初期条件と時刻から直交座標系での状態ベクトルを計算<br>
    •概要<br>
    時刻と離心近点離角等の変数を結びつける超越方程式(ケプラー方程式やバーカー方程式)を解き,<br>
    初期条件と与えられた時刻から直交座標系での状態ベクトルを計算する.<br>
    アルゴリズムは3種類あり，<br>
    a) 状態ベクトルを焦点中心座標系の単位ベクトルで表現するもの<br>
    b) ラグランジュの係数を用いて表現するもの<br>
    c) 普遍変数を導入した上でラグランジュの係数を用いて表現するもの<br>
    aのアルゴリズムについて，Curtis(2010)では地球の扁平率を考慮した方法が記されているため，"practice4ac"として作成した．
