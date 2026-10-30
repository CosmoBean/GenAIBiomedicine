% cisXXX.cls -- minimal homework class with the same interface as the
% UPenn CS HW Overleaf template (\HWauthor, \HWno, \HWcourse, \HWpartner,
% \HWproblem, \HWsubproblem). If you have the original cisXXX.cls from
% Overleaf, it can replace this file without changing report.tex.
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cisXXX}[2026/10/29 homework class]
\LoadClass[11pt]{article}
\RequirePackage[margin=1in,headheight=14pt]{geometry}
\RequirePackage{fancyhdr}
\RequirePackage{xcolor}

\newcommand{\@hwauthor}{}\newcommand{\@hwemail}{}
\newcommand{\@hwno}{}\newcommand{\@hwcourse}{}
\def\@hwpartners{}
\newcommand{\HWauthor}[2]{\renewcommand{\@hwauthor}{#1}\renewcommand{\@hwemail}{#2}}
\newcommand{\HWno}[1]{\renewcommand{\@hwno}{#1}}
\newcommand{\HWcourse}[1]{\renewcommand{\@hwcourse}{#1}}
\newcommand{\HWpartner}[1]{\g@addto@macro\@hwpartners{\par #1}}

\pagestyle{fancy}
\fancyhf{}
\lhead{\@hwauthor}
\chead{\@hwcourse}
\rhead{Homework \@hwno}
\cfoot{\thepage}

\renewcommand{\maketitle}{%
  \thispagestyle{fancy}%
  \begin{center}
    {\Large\bfseries \@hwcourse\ \,--\,\ Homework \@hwno\par}\vspace{4pt}
    {\large \@hwauthor\par}
    {\ttfamily \@hwemail\par}
    \ifx\@hwpartners\@empty\else{\small Collaborators:\@hwpartners\par}\fi
  \end{center}\vspace{2pt}\hrule\vspace{8pt}}

\newcounter{HWproblem}
\newcounter{HWsubproblem}[HWproblem]
\renewcommand{\theHWsubproblem}{\alph{HWsubproblem}}
\newcommand{\HWproblem}[1][]{%
  \par\bigskip\stepcounter{HWproblem}%
  \noindent{\large\bfseries Problem \theHWproblem\if\relax\detokenize{#1}\relax\else\ --\ #1\fi}%
  \par\nobreak\vspace{2pt}\hrule\nobreak\medskip}
\newcommand{\HWsubproblem}[1][]{%
  \par\medskip\refstepcounter{HWsubproblem}%
  \noindent\textbf{(\theHWsubproblem)}\if\relax\detokenize{#1}\relax\else\ \textbf{#1}\fi\quad}
