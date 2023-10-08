function debugThis(t, tk) {
    if (t >= tk) {
        debugger
    }
}
function printLog(t, el) {
    console.log(`t = ${t.toFixed(5)} m = ${el[0].toFixed(4)} vx = ${el[1].toFixed(4)} vy = ${el[2].toFixed(4)} x = ${el[3].toFixed(4)} y = ${el[4].toFixed(4)} v = ${el[5].toFixed(4)} r = ${el[6].toFixed(4)} `);
}
export {
    debugThis,
    printLog
}