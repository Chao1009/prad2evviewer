#include <QApplication>
#include <QWebEngineView>
#include <QCommandLineParser>
#include <QUrl>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("prad2_monitor");
    app.setApplicationVersion("1.0.0");

    QCommandLineParser parser;
    parser.setApplicationDescription("Lightweight Qt WebEngine client for PRad2 event monitor");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("url", "URL to load (default: http://localhost:8080)",
                                "[url]");

    parser.process(app);

    QString url = "http://localhost:8080";
    const QStringList args = parser.positionalArguments();
    if (!args.isEmpty()) {
        url = args.first();
    }

    QWebEngineView view;
    view.setWindowTitle("PRad2 Monitor");
    view.resize(1280, 800);
    view.load(QUrl(url));
    view.show();

    return app.exec();
}
