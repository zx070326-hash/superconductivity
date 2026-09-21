---
name: deliver
description: 收尾，准备交给用户。做动作、不印清单：跑门禁、列档位、看 diff。
---

交付物：`$ARGUMENTS`。执行并把输出贴进交付：

1. `git status --short` 与 `git diff --stat`：只动了该动的文件？分支 UNMERGED？未 commit / push？
2. `sh tools/check_cited_scripts_tracked.sh`：被笔记点名的脚本都已 `git add -f`？
3. 新建的每个 `theory/notes/**/*.md` 与 `theory/derivations/**/*.md`：在对应 INDEX.md 里有行？取代旧笔记的，旧行写 `superseded by <文件名>`。
4. 逐条列本次结论 + 档位（analytic / rational-exact 含格数与有理数个数、不得简称 analytic / certificate / numerical-evidence）。越级的降回去。
5. 承重项：哪几个是我从零重算的，哪几个是转引的。转引的写明来源，不写成已核。
6. 若含订正：新陈述的第二段（凭什么对）写了吗？
7. 若含收紧：预算与消费点写了吗？∞ 的那些标为「只买余量」。
8. 首句写没做完/失败/跳过的部分；没有就写「全部完成」。

本条可被覆盖，写明理由即可。
