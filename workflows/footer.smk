## Queue-agnostic fallback resources for all rules.
## These defaults apply when a rule has no explicit resource set.
## They keep local/non-Slurm runs from lacking baseline resource values.
_fallback_rule_resources = {
	"mem_mb": 20000,
	"runtime": 480,
	"nodes": 1,
}

_user_rule_defaults = config.get("DEFAULT_RULE_RESOURCES", {})
if isinstance(_user_rule_defaults, dict):
	_fallback_rule_resources.update(_user_rule_defaults)


def _has_resource(rule_obj, resource_name):
	try:
		return resource_name in rule_obj.resources.keys()
	except Exception:
		return hasattr(rule_obj.resources, resource_name)


def _set_resource(rule_obj, resource_name, value):
	try:
		rule_obj.resources[resource_name] = value
	except Exception:
		setattr(rule_obj.resources, resource_name, value)


for _rule in workflow.rules:
	for _res_name, _res_value in _fallback_rule_resources.items():
		if not _has_resource(_rule, _res_name):
			_set_resource(_rule, _res_name, _res_value)

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
