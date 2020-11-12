"""
Perhaps, although functionality can be here, move authentication/authorization
to NGINX so it can filter ALL requests. An example seemingly interesting project:
https://github.com/mbreese/subauth

"""
from pyparsing import ParserElement, Literal, oneOf, CaselessKeyword, Word, alphas, alphanums, QuotedString, \
    delimitedList, Group, operatorPrecedence, Or, opAssoc, Dict

import biobarcoding
from biobarcoding.db_models.sysadmin import Identity, Role, Organization, RoleIdentity, GroupIdentity, \
    OrganizationIdentity

ParserElement.enablePackrat()
lparen, rparen, lbracket, rbracket, lcurly, rcurly, dot, equals, hash = map(Literal, "()[]{}.=#")
double_quote = Literal('"')
single_quote = Literal("'")
quote = oneOf('" ''')  # Double quotes, single quotes
comparisonop = oneOf("= == != <>")
andop = CaselessKeyword("AND")
orop = CaselessKeyword("OR")
notop = CaselessKeyword("NOT")
inop = CaselessKeyword("IN")
simple_ident = Word(alphas, alphanums+"_")  # Start in letter, then "_" + letters + numbers
in_expression = Group(simple_ident + inop.suppress() + lparen.suppress() + delimitedList(Or([QuotedString('"'), QuotedString("'")]), ",") + rparen.suppress()).\
    setParseAction(lambda t: {'type': 'in',
                              'attribute': t[0][0],
                              'values': t[0][1:]}
                   )

authr_expression = operatorPrecedence(Or([in_expression, simple_ident, QuotedString('"'), QuotedString("'")]),
                                      [
                                         (comparisonop, 2, opAssoc.LEFT, lambda _s, l, t: {
                                             'type': 'comparison',
                                             'terms': t.asList()[0][0::2],
                                             'ops': t.asList()[0][1::2]
                                         }),
                                         (notop, 1, opAssoc.RIGHT, lambda _s, l, t: {
                                             'type': 'not',
                                             'terms': [t.asList()[0][1]],
                                             'ops': [t.asList()[0][0]]
                                         }),
                                         (andop, 2, opAssoc.LEFT, lambda _s, l, t: {
                                             'type': 'and',
                                             'terms': t.asList()[0][0::2],
                                             'ops': t.asList()[0][1::2]
                                         }),
                                         (orop, 2, opAssoc.LEFT, lambda _s, l, t: {
                                             'type': 'or',
                                             'terms': t.asList()[0][0::2],
                                             'ops': t.asList()[0][1::2]
                                         }),
                                      ],
                                      lpar=lparen.suppress(),
                                      rpar=rparen.suppress()
                                      )


def string_to_ast(rule: ParserElement, input_: str) -> Dict:
    """
    Convert the input string "input_" into an AST, according to "rule"

    :param rule:
    :param input_:
    :return: a dictionary conforming the AST (the format changes from rule to rule)
    """

    def clean_str(us):
        # "En dash"                 character is replaced by minus (-)
        # "Left/Right double quote" character is replaced by double quote (")
        # "Left/Right single quote" character is replaced by single quote (')
        # "€"                       character is replaced by "eur"
        # "$"                       character is replaced by "usd"
        return us.replace(u'\u2013', '-'). \
            replace(u'\u201d', '"').replace(u'\u201c', '"'). \
            replace(u'\u2018', "'").replace(u'\u2019', "'"). \
            replace('€', 'eur'). \
            replace('$', 'usd')

    res = rule.parseString(clean_str(input_), parseAll=True)
    res = res.asList()[0]
    while isinstance(res, list):
        res = res[0]
    return res


# Comparison operators
op_map = {
        "==": lambda a, b: a == b,
        "=": lambda a, b: a == b,
        "!=": lambda a, b: a != b,
        "<>": lambda a, b: a != b,
        }


def ast_evaluator(ast: Dict, ident: Identity):
    if "type" in ast:
        if ast["type"] == "in":
            if ast["attribute"] in ("role", "group", "organization", "identity"):
                if ast["attribute"] in ("role", "group", "organization"):
                    lst_names = [getattr(o, ast['attribute']).name for o in getattr(ident, f"{ast['attribute']}s")]
                else:
                    lst_names = [ident.name]
                found = False
                for name in lst_names:
                    if name in ast["values"]:
                        found = True
                        break
                return found
            else:
                raise Exception("Attribute not supported")
        elif ast["type"] in ("not", "and", "or", "comparison"):
            current = None
            for i, e in enumerate(ast["terms"]):
                following = ast_evaluator(e, ident)
                if i == 0:
                    current = following
                    op = ast["ops"][i].lower()
                    if op == "not":
                        current = not bool(following)
                else:
                    op = ast["ops"][i-1].lower()
                    if op == "and":
                        current = current and following
                    elif op == "or":
                        current = current or following
                    else:  # Comparators
                        fn = op_map[op]
                        current = fn(current, following)
            return current
        else:
            raise Exception(f"Not supported type {ast['type']}")


if __name__ == '__main__':
    ident = Identity()
    role = Role()
    group = biobarcoding.db_models.sysadmin.Group()
    organization = Organization()
    ident.name = "rnebot"
    role.name = "admin"
    group.name = "devs"
    organization.name = "itc"
    role_id = RoleIdentity()
    role_id.identity = ident
    role_id.role = role
    group_id = GroupIdentity()
    group_id.identity = ident
    group_id.group = group
    org_id = OrganizationIdentity()
    org_id.identity = ident
    org_id.organization = organization

    in_examples = ["role in ('admin', 'advanced_user')",
                   "group in ('nextgendem', 'students')"]

    full_examples = ["(role in ('admin', 'advanced_user') and group in ('nextgendem')) or not role in ('people') and not group in ('students')"]

    for e in in_examples:
        # try:
            ast = string_to_ast(in_expression, e)
            print(ast)
            res = ast_evaluator(ast, ident)
            print(res)
        # except:
        #     print("Incorrect")

    for e in full_examples:
        # try:
            print(e)
            ast = string_to_ast(authr_expression, e)
            print(ast)
            res = ast_evaluator(ast, ident)
            print(res)
        # except:
        #     print("Incorrect")