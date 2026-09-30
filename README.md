# RockSolver V0–V6

RockSolver 是独立于 Vue、PyQt、OpenSees 和 FLAC3D 的 C++20 岩体计算核心。同一份 `rocksolver_core` 被 Native CLI 与 Emscripten WASM adapter 共用；Web 只生成 `RockModel`、调度求解和渲染 `RockResult`。

当前纵切覆盖：

- V0：CMake、RockModel/RockResult、Native/WASM adapter、原生测试；
- V1：规则网格拆分 Tet4、小应变线弹性、重力、边界和显式动态松弛；
- V2：与 `cmodelsmrm007` 同更新次序的 SMRM 材料点：完整 6×6 各向异性柔度、三组结构面、应力相关闭合 `k`、剪切退化 `h`、增量应力更新和完整诊断量；
- V3：参数化边坡与强度折减阶段；
- V4：最多三组结构面及开闭、滑移、拉/剪破坏状态；
- V5：隧道初始应力、分段开挖失活、围岩松弛和等效支护弹簧；
- V6：版本化 JSON 协议、63 项 FLAC 兼容诊断输出、资源估算、本地 HTTP 服务、跨平台构建入口。

此外，Native CLI、本地 HTTP 服务和浏览器 WASM 都支持导入 Abaqus `C3D4` 网格。网页端可以直接打开原 FLAC3D 建模完整包 ZIP：前端提取 `*.inp` 与命令流参数，浏览器 Web Worker 默认使用原始节点和四面体单元计算，不再重新生成规则网格，也不依赖本机服务。

## 主结果格式

新计算结果以 `.rockresult` 作为主文件。它是可版本化的压缩结果容器，包含：

- `manifest.json`：格式版本、求解器版本、节点/单元数量、工程关联和后处理算法；
- `model.json`：材料、结构面、阶段和求解参数，不重复保存原始 INP 文本；
- `result.json`：网格、位移、应力、应变、本构状态和 JavaScript SMRM 后处理结果；
- `README.txt`：格式用途说明。

参数场直接读取该结果包，不经过 CSV 中转。CSV 只由网页按需导出，用于 Excel、FLAC3D 对比和第三方交换。结果包清单规范见 `schemas/rock-result-package.schema.json`。

材料点算法已经按原 `cmodelsmrm007` 的应力符号、Voigt 顺序、结构面方向、柔度矩阵和旧状态增量更新次序实现，并保留路径历史。RockSolver 的全局方程仍是 Tet4 有限元显式动态松弛，不是 FLAC3D 的有限差分/有限体积 Zone 求解器；因此不能把整场结果描述为逐点数值等价。当前支护也不是独立 Bolt/Shell/Lining 单元。工程使用前仍应按 `docs/VALIDATION_PLAN.md` 逐级与 FLAC3D 结果对照。

## Native

```bash
./scripts/build_rocksolver_macos.sh
PYTHONPATH=src .venv/bin/python -m geostudio.rocksolver.service
```

直接检查或计算一个 INP：

```bash
runtime/rocksolver/macos/arm64/rocksolver_native --estimate --inp /path/to/your_model.inp
runtime/rocksolver/macos/arm64/rocksolver_native --scenario slope --material smrm --inp /path/to/your_model.inp
```

目前导入器只接受三维四节点四面体 `C3D4`。它保留原节点和单元拓扑，并按当前两个完整包的命令建立底面 XYZ、X 外侧法向和 Y 外侧法向约束。原命令流中超出已识别子集的复杂分组、接触、渗流、动力、开挖和实体支护命令不会被臆测执行，需扩展协议并增加回归算例。

SMRM 结构面的抗拉强度为零时，任意正法向应力都会立即判定为张拉破坏。为避免 RR 输出为无穷大并破坏 JSON，严重度使用完整岩石抗拉强度作有限归一化；这只影响诊断数值的尺度，不会给零抗拉结构面增加承载力。

Windows 使用 `scripts/build_rocksolver_windows.ps1`。服务仅默认监听 `127.0.0.1:8766`。

## WASM

激活 Emscripten SDK 后运行：

```bash
./scripts/build_rocksolver_wasm.sh
```

输出 `rocksolver/build/wasm-release/rocksolver_wasm.js` 和 `.wasm`。传入 Vue 的 `public/rocksolver` 作为脚本参数即可安装 Web 资产。Native 与 WASM 都链接同一个 `rocksolver_core` 静态库，WASM adapter 支持参数化模型和完整 INP 模型。
