# YC38 publication and original local history

10 October 2026. Published with the owner's explicit approval.

Direct Git push was unavailable because the runtime had no GitHub CLI
credential. The connected GitHub app publishes the complete YC38 working
tree relative to the existing compact-observer branch, plus the original
Git history in `yc38-local-history.bundle`. This is a snapshot publication,
not a claim that the remote snapshot commit has the old local commit ID.

`YC38_PUBLICATION.json` records the original head/tree, 50 previously
unpublished commits, bundle prerequisites and SHA-256. The bundle preserves
the original commit identities used in frozen source manifests. To recover
them in a clone containing the prerequisites:

```sh
git bundle verify research-history/yc38-local-history.bundle
git fetch research-history/yc38-local-history.bundle refs/heads/research/yc38-certificate-audit-2026-10-10:refs/heads/archive/yc38-original
```

Existing source and certificate files are copied byte-for-byte from the
local YC38 tree. No historical certification result is newly replayed or
promoted by this publication. Further theoretical stages remain separate.
