"""Custom Pygments styles for the multicolorfits docs code blocks.

Same denim light/dark machinery as skyplothelper (gedit schemes); site chrome
uses mcf accents in ``_static/custom.css``.
"""

from pygments.style import Style
from pygments.token import (
    Comment,
    Error,
    Generic,
    Keyword,
    Literal,
    Name,
    Number,
    Operator,
    Punctuation,
    String,
    Token,
)


class DenimDarkStyle(Style):
    """Dark-mode code style (gedit denimdark)."""

    name = "denimdark"
    background_color = "#211f1d"  # warm charcoal, matches the dark surface
    highlight_color = "#3f3f3f"

    styles = {
        Token: "#d9d5c5",              # cream text
        Comment: "italic #7e736a",     # warm grey
        Comment.Preproc: "#7e736a",
        Comment.Special: "italic #7e736a",
        Keyword: "#597F88",            # denim
        Keyword.Constant: "#b54218",   # rust (True/False/None)
        Keyword.Type: "bold #597F88",
        Operator: "#49688D",
        Operator.Word: "#49688D",
        Punctuation: "#d9d5c5",
        Name: "#d9d5c5",
        Name.Builtin: "#49688D",
        Name.Builtin.Pseudo: "#b54218",
        Name.Function: "#455466",
        Name.Function.Magic: "#455466",
        Name.Class: "bold #455466",
        Name.Namespace: "#49688D",
        Name.Exception: "#b54218",
        Name.Decorator: "#9b703f",     # copper
        Name.Constant: "#9b703f",
        Name.Variable: "#d9d5c5",
        Name.Attribute: "#d9d5c5",
        Name.Tag: "#49688D",
        String: "#8f9d6a",             # sage
        String.Doc: "italic #8f9d6a",
        String.Escape: "#9b703f",
        String.Interpol: "#9b703f",
        Number: "#9b703f",             # copper
        Literal: "#9b703f",
        Generic.Heading: "bold #c5b777",    # gold
        Generic.Subheading: "bold #c5b777",
        Generic.Emph: "italic",
        Generic.Strong: "bold",
        Generic.Deleted: "#b54218",
        Generic.Inserted: "#8f9d6a",
        Generic.Output: "#7e736a",
        Generic.Prompt: "#49688D",
        Error: "#d9d5c5 bg:#551212",
    }


class DenimLightStyle(Style):
    """Light-mode code style (gedit denimlight)."""

    name = "denimlight"
    background_color = "#efe7d6"  # warm parchment, matches the light header band
    highlight_color = "#e3dcc9"

    styles = {
        Token: "#383838",             # near-black ink
        Comment: "italic #8f8f8f",    # grey
        Comment.Preproc: "#8f8f8f",
        Comment.Special: "italic #8f8f8f",
        Keyword: "#2e5266",           # denim
        Keyword.Constant: "#b54218",  # rust
        Keyword.Type: "bold #2e5266",
        Operator: "#2e5266",
        Operator.Word: "#2e5266",
        Punctuation: "#383838",
        Name: "#383838",
        Name.Builtin: "#2e5266",
        Name.Builtin.Pseudo: "#b54218",
        Name.Function: "#2e5266",
        Name.Function.Magic: "#2e5266",
        Name.Class: "bold #2e5266",
        Name.Namespace: "#2e5266",
        Name.Exception: "#b54218",
        Name.Decorator: "#7d5a33",    # brown
        Name.Constant: "#7d5a33",
        Name.Variable: "#383838",
        Name.Attribute: "#383838",
        Name.Tag: "#2e5266",
        String: "#759568",            # green
        String.Doc: "italic #759568",
        String.Escape: "#7d5a33",
        String.Interpol: "#7d5a33",
        Number: "#7d5a33",            # brown
        Literal: "#7d5a33",
        Generic.Heading: "bold #2e5266",
        Generic.Subheading: "bold #2e5266",
        Generic.Emph: "italic",
        Generic.Strong: "bold",
        Generic.Deleted: "#b54218",
        Generic.Inserted: "#759568",
        Generic.Output: "#8f8f8f",
        Generic.Prompt: "#2e5266",
        Error: "#383838 bg:#f0d2d2",
    }
