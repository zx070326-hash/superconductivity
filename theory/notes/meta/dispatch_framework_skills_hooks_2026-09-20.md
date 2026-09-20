# 派活框架：skill + hook 机制（2026-09-20 本机实测）

这套东西解决**多天或多轮协作中，纪律预期不在脑子里**的问题。事故清单：派活未算预算（三起）、订正时新陈述被免检放行（W1 失效热区，一轮四起）、未索引笔记（一天三篇 + deficit-real-solve 两篇其一已被自己的后一次交付取代）、agent 交付三天前就在仓库里（重复交付）、引用脚本门禁假 PASS 一整天（不递归 + 正则无连字符 ⟹ 09-18 重组后扫不到任何笔记）。

## 本机实测机制（2026-09-20，三只 haiku 金丝雀 + 一个探针 skill）

- **PreToolUse(Agent)** 的 `additionalContext` 只有 lead 看得到，子 agent 看不到。
- `updatedInput.prompt` 真的改写子 agent 收到的 prompt（子 agent 原样抄出被前置的那行）。
- 两者可在同一份 `hookSpecificOutput` 里并存。
- 子 agent 同时加载 `CLAUDE.md` 与 `CLAUDE.local.md` ⟹ 派活 prompt 里不再抄执行方纪律。
- 项目 skill 在 `.claude/skills/<name>/SKILL.md`，`description` 逐字进入模型可用列表（新建后下一轮刷新），正文只在调用时进上下文，`$ARGUMENTS` 原样替换。
- Windows 下 hook 用 Python + `sys.stdin.buffer` + `ensure_ascii`。测试 payload 必须用 `json.dumps` 生成，`printf` 会折掉反斜杠。

## hook_prior_art.py（PreToolUse: Agent）

- 触发：`Agent` 工具被调用时。
- 行为：抽取 `description + prompt` 文本 ⟹ 管道喂给 `prior_art.sh` ⟹ 追加到子 agent 的 prompt 末尾（`updatedInput.prompt` 改写）、同时给 lead 看一份（`additionalContext`）。
- 不阻断、不判决。任何异常 ⟹ 静默放行（exit 0、无输出）。
- `prior_art.sh` 以 `tools/vocab.tsv` 为受控词表，对承重量做固定串匹配（`grep -F`），按「同时命中几个承重量」排序，分流 live / 非 live，标行内撤回。输出必须短（≤ 20 行），进主会话上下文。
- 输出追加块显式标明来源（「由 hook 生成，不是派活方的判断」），执行方知道那不是派活方写的。
- 先验检索不是判决，重算作为独立二次验证是正当的；本条只是确保「已存在」不被忽略。

## hook_note_index.py（PostToolUse: Write|Edit）

- 触发：`theory/notes/**/*.md` 或 `theory/derivations/**/*.md` 被 Write/Edit，且文件名不在对应 `INDEX.md` 里（检测链接模式 `[<文件名>](`）。
- 行为：只提醒（`additionalContext` 给当前会话），不阻断、不改文件。
- 内容：「路由规则要求落盘同一次编辑里补索引行（含状态列；若取代旧笔记，旧行写 `superseded by <文件名>`）。」
- 任何异常 ⟹ 静默放行。

## 九个 skill（什么时候选它 / 行数）

| skill | 选它的时刻（与最近邻的区分，逐字取自 description） | 行数 |
|---|---|---|
| `dispatch-tighten` | 派活，且交付物会被某个判决直接消费（收紧/认证/复核一个承重量）。产出一份 6 格 prompt。填不出预算或消费点 ⟹ 不派。 | 21 |
| `dispatch-explore` | 派活，但预算填不出来（侦察文献、扫可能性、找反例形状）。选它即声明「探索型，预算门不适用」。产出一份 3 格 prompt。 | 17 |
| `receive` | 收活：子 agent 交付回来了，我在首次接收与分诊、尚未下判断。派活与收活隔了几小时，派活时的判据已不在上下文里。要落笔说「它错了」时转 verify-correction。 | 17 |
| `verify-correction` | 我正要落笔说「这个错了」——对象是别人或我自己的旧陈述（首次接收分诊在 receive）。W1 失效热区：注意力全在旧错，新陈述被免检放行。 | 23 |
| `tighten-bound` | 我自己在纸笔/符号推导里放缩或收紧一个量，产出是不等式。不是派给别人（dispatch-*），不是写区间算术脚本（write-certificate）。 | 16 |
| `design-test` | 我自己在设计数值检验（抽样/扫描/对照），产出是通过率或反例，不要求声音性（float64 即可）。与 tighten-bound 的区别：造证据，不证不等式；与 write-certificate 的区别：产物不是区间包围。 | 16 |
| `write-certificate` | 在写或改区间算术/证书脚本（iv/mpmath），产出是声音的上/下界。与 design-test 的区别：这里的产物要求声音性，验收门是投毒必须报警。 | 16 |
| `read-source` | 承重字符（公式、表格数、约定）来自外部 PDF/网页。纯机械清单。 | 17 |
| `deliver` | 收尾，准备交给用户。做动作、不印清单：跑门禁、列档位、看 diff。 | 17 |

说明：
- 一个 skill 对应一个「纪律预期不在脑子里」的时刻，不对应一条纪律。
- 机制能自动产出的不做 skill（做 hook）。
- ≤20 行、只说查什么不说得出什么、末行可覆盖。
- 多方向派活不设单独 skill（= 同一 skill × N + 一行共享检查）、收活单独设。

## 设计规则：刻意不建的

- `experiments/` 文件夹 — 状态已在 `INDEX.md` 状态列，不另立目录分层。
- 根 `CLAUDE.md` 新增条目 — 本文件只管「怎么做」，演进中的状态在 auto-memory。
- `CLAIMS.tsv` 账本 — 复用 `INDEX.md` 已有的状态列，避免两处真源漂移 + 台账腐烂。
- 12 份子目录宪法 — 只建两份（`theory/notes/CLAUDE.md`、`theory/derivations/CLAUDE.md`），内容是从两个 INDEX 头部**搬**来的既有约定，不立新法；INDEX 头部只留指针。
- `dispatch_template.md` — 已删，正文并入 skill。模板是静态的，纪律是从事故积累的；前者腐烂快、后者活在 skill 里。

## 对抗审查（2026-09-20，sonnet reviewer）

采纳 8 条：tighten-bound/design-test/receive/verify-correction 的 description 互相点名最近邻；rational-exact 三处补「不得简称 analytic」；receive 新增「界收紧/判决翻转」一行（数字错 vs 证书松、成本单位、重解下游）；design-test 自由度补「找生成式」；hook 失效 fallback；引文未见不得承重。
拒绝 5 条：dispatch-tighten 容差写法（已在第 13 行，审查漏看）；原子型概率公式（skill 只说查什么）；design-test 加宽度断言（它不产区间）；标级三处重复（有意）；「配对生成检验」改泛化（有意保留具体案例）。

## 状态

全部 UNCOMMITTED 在 `p4-conj2`，等用户放行。包括：
- `.claude/settings.json` 两个 hook 配置
- `.claude/skills/*/SKILL.md` 九个 skill
- `tools/hook_prior_art.py` / `tools/hook_note_index.py` / `tools/prior_art.sh` / `tools/vocab.tsv`
- 本笔记
